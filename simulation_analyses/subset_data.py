from Bio import SeqIO

input_file = "/project/mpennell_978/kavoss/real_data_cf/cattle/14363/1/1/day28_z/clean.fasta"
output_file = "/project/mpennell_978/kavoss/real_data_cf/cattle/14363/1/1/day28_z_subset/clean.fasta"

# Load sequences from the input fasta file
sequences = list(SeqIO.parse(input_file, "fasta"))

# Calculate 20% of the total sequences
total_sequences = len(sequences)
subset_count = int(0.2 * total_sequences)

print(f"Total sequences: {total_sequences}")
print(f"Extracting first {subset_count} sequences...")

# Write the first 20% of sequences to the output file
SeqIO.write(sequences[:subset_count], output_file, "fasta")

print(f"Subset file created: {output_file}")
