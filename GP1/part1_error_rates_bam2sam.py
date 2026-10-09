# BAM to SAM
import argparse
import subprocess

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
		sequence = fields[9] 	# read bases column in SAM file

		# md tags
		md = None
		for tag in fields[11:]:
			if tag.startswith("MD:Z:"):
				md = tag[5:]	# removes MD:Z:
				break
		
		if md is None:
			raise ValueError(f"Read {fields[0]} is missing an MD tag")
		
'''		# sample test first three primary mapped aligns
		if primary_count <= 3:
			print(f"Read: {fields[0]}")
			print(f"SEQ (first 30 bases): {sequence[:30]}")
			print(f"CIGAR: {cigar}")
			print(f"MD: {md}")
			print()

print(f"Primary mapped alignments: {primary_count}")
'''

