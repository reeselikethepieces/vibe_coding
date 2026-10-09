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