from Bio import AlignIO
import sys

# Specify the input and output file names
input_fasta = sys.argv[1]
output_phylip = sys.argv[2]

# Convert the FASTA alignment to PHYLIP format
AlignIO.convert(input_fasta, "fasta", output_phylip, "phylip")

print(f"Successfully converted {input_fasta} to PHYLIP format.")
