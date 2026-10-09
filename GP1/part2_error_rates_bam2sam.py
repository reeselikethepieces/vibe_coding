# BAM to SAM
import argparse
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
			
parser = argparse.ArgumentParser(description="Convert a BAM file to SAM")
parser.add_argument("bam", help="Input BAM file")
parser.add_argument("-o", "--output", required=True, help="Output SAM file")
args = parser.parse_args()

subprocess.run(
	["samtools", "view", "-h", "-o", args.output, args.bam],
	check=True,
)

# read SAM files and skip unmapped, secondary, and supplementary alignments
primary_count = 0

with open(args.output, "r") as sam_file:
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

# print(f"Primary mapped alignments: {primary_count}")
# can check above value by using 'samtools view -c -F 2308 file.ext'

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
		
		if primary_count <= 3:
			print(f"Read: {fields[0]}")
			print(f"Evaluated positions: {len(base_pairs)}")
			print(f"First 10 pairs (reference, read): {base_pairs[:10]}")
			print()

print(f"Primary mapped alignm,ents: {primary_count}")
'''		# sample test first three primary mapped aligns
		if primary_count <= 3:
			print(f"Read: {fields[0]}")
			print(f"SEQ (first 30 bases): {sequence[:30]}")
			print(f"CIGAR: {cigar}")
			print(f"MD: {md}")
			print()

print(f"Primary mapped alignments: {primary_count}")
# next is to add a function that returns pairs and to import re
'''