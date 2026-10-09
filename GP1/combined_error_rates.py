import argparse
import csv
from pathlib import Path
import re
import subprocess

def get_base_pairs(seq, cigar, md):
	# return (ref base, read base) pairs
	seq = seq.upper()
	
	if seq == "*" or "=" in seq:
		raise ValueError("Explicit read bases are required")
	if not re.fullmatch(r"(?:[1-9][0-9]*[MIDNSHP=X])+", cigar):
		raise ValueError("Invalid CIGAR")
	if not re.fullmatch(r"[0-9]+(?:(?:[A-Z]|\^[A-Z]+)[0-9]+)*", md):
		raise ValueError('Invalid MD tag')
	
	# build alignment columns from MD
	columns = []
	read_index = 0
	
	for length, operation in re.findall(r"([0-9]+)([MIDNSHP=X])", cigar):
		length = int(length)
		
		# below 'in' is shorthand for - if operation = 'M' or ... 
		if operation in "M=X":
			bases = seq[read_index:read_index + length]
			columns.extend((base, operation) for base in bases)
			read_index += length
		# below 'in' is shorthand for - elif operation == "I" or operation == "S":
		elif operation in "IS":
			read_index += length
		elif operation == "D":
			columns.extend((None, "D") for _ in range (length))
		
		# H, N, and P in CIGARs consume no read bases and are absent from MD

	if read_index != len(seq):
		raise ValueError('CIGAR length does not match the read seq')
	
	# walk through MD and alignment columns together
	pairs = []
	column_index = 0
	
	for token in re.findall(r"[0-9]+|\^[A-Z]+|[A-Z]", md):
		if token.isdigit():
			event, length = "match", int(token)
		elif token.startswith("^"):
			event, length = "deletion", len(token) - 1
		else:
			event, length = "mismatch", 1
		
		if column_index + length > len(columns):
			raise ValueError("MD extends beyond CIGAR")
		
		for _ in range(length):
			read_base, operation = columns[column_index]
			column_index += 1
			
			if event == "deletion":
				if operation != "D":
					raise ValueError("MD deletion conflicts with CIGAR")
				continue
			
			if event == "match":
				if operation not in "M=":
					raise ValueError("MD mismatch conflicts with CIGAR")
				reference_base = read_base
			
			else:
				if operation not in "MX" or read_base == token:
					raise ValueError("MD mismatch conflicts with CIGAR or SEQ")
				reference_base = token
			
			# keep only positions with definite A/C/G/T bases
			if reference_base in "ACGT" and read_base in "ACGT":
				pairs.append((reference_base, read_base))

	if column_index != len(columns):
		raise ValueError("MD ends before CIGAR")
	
	return pairs
			
parser = argparse.ArgumentParser(description="Pool 12 substitution rates across input BAM files")
parser.add_argument("bam", nargs="+", help="one or more BAM files")
parser.add_argument("-o", "--output", help="SAM filename when analyzing a single BAM")
parser.add_argument("--sam-dir",default=".", help="Directory for converted SAM files (default: current directory)")
parser.add_argument("--tsv", required=True, help="output pooled sugbsitution-rate TSV file")
args = parser.parse_args()

# -o to name SAM output file and ensures one BAM files
if args.output and len(args.bam) != 1:
	parser.error("-o names a single SAM; for multiple BAMs use --sam-dir instead")

# prepare and check paths for multiple BAM files

# build list of full BAM paths, expanduser() expands ~; resolve() makes each path absolute
bam_paths = [Path(name).expanduser().resolve() for name in args.bam]

if len(set(bam_paths)) != len(bam_paths):
	parser.error("The same BAM was supplied more than once")

for bam in bam_paths:
	if not bam.is_file():
		parser.error(f"Input BAM does not exist: {bam}")

# determine where converted SAM files will go
sam_dir   = Path(args.sam_dir).expanduser().resolve()
sam_paths = ([Path(args.output).expanduser().resolve()]
			if args.output
			else [sam_dir / (bam.stem + ".sam") for bam in bam_paths])
tsv_path  = Path(args.tsv).expanduser().resolve()
all_paths = bam_paths + sam_paths + [tsv_path]

if len(set(all_paths)) != len(all_paths):
	parser.error("All input and output paths must be different; BAMs also need distinct names")
sam_dir.mkdir(parents=True, exist_ok=True)

# read SAM files and skip unmapped, secondary, and supplementary alignments
primary_count = 0
BASES = "ACGT"

# list comprehension
substitution_types = [
	(ref, read)			# tuple
	for ref in BASES
	for read in BASES
	if ref != read
]

# dictionaries
reference_counts = {base:0 for base in BASES}
substitution_counts = {pair: 0 for pair in substitution_types}

# shared counters above this loop accumulate counts from every BAM
for bam_path, sam_path in zip(bam_paths, sam_paths):
	subprocess.run(
		["samtools", "view", "-h", "-o", str(sam_path), str(bam_path)],
		check=True
	)
	file_primary_count = 0
	with open(sam_path, "r") as sam_file:
		for line in sam_file:
			# skips header lines
			if line.startswith("@"):
				continue
		
			fields = line.rstrip("\n").split("\t")
			flag = int(fields[1])
		
			# skip unmapped, secondary, or supp aligns
			if flag & (0x4 | 0x100 | 0x800):
				continue
		
			primary_count += 1
			file_primary_count += 1
			cigar = fields[5] 		# alignment structure column in SAM file
			seq = fields[9] 	# read bases column in SAM file

		# md tags
			md = None
			for tag in fields[11:]:
				if tag.startswith("MD:Z:"):
					md = tag[5:]	# removes MD:Z:
					break
		
			if md is None:
				raise ValueError(f"Read {fields[0]} is missing an MD tag")

			base_pairs = get_base_pairs(seq, cigar, md)
		
			for reference_base, read_base in base_pairs:
				reference_counts[reference_base] += 1
			
				if reference_base != read_base:
					substitution_counts[(reference_base, read_base)] += 1
		# above is crucial because we consider all evaluated pairs as part of the total
	print(f"{bam_path.name}: {file_primary_count} primary mapped alignments")

print(f"Primary mapped alignments: {primary_count}")

# these are pooled counts divided by pooled ref totals not average of rates
rates = []

for ref, read in substitution_types:
	total = reference_counts[ref]
	if total > 0:
		rate = substitution_counts[(ref, read)] / total
		rates.append(f"{rate:.12g}")
	else:
		rates.append("NA")

rate_lookup = dict(zip(substitution_types, rates))

with open(args.tsv, "w", newline="") as tsv_file:
	writer = csv.writer(tsv_file, delimiter="\t", lineterminator="\n")
	writer.writerow(["reference/read"] + list(BASES))
	
	for ref in BASES:
		row = [
			"-" if ref == read else rate_lookup[(ref,read)]
			for read in BASES
		] 
		writer.writerow([ref] + row)

print(f"Substitution matrix written to: {args.tsv}")