import pandas as pd
import os
import sys

path = sys.argv[1]
output = sys.argv[2]
# Read the CSV file into a pandas DataFrame
df = pd.read_csv(path)

# Group DataFrame by TimePoint
grouped = df.groupby('TimePoint')

# Output directory for FASTA files
output_dir = output
os.makedirs(output_dir, exist_ok=True)

# Iterate over each TimePoint group
for timepoint, group in grouped:
    # Create a FASTA file for the current TimePoint
    fasta_filename = os.path.join(output_dir, f'{timepoint}.fasta')
    with open(fasta_filename, 'w') as fasta_file:
        # Iterate over rows in the current group
        for index, row in group.iterrows():
            # Write sequence header and sequence to the FASTA file
            fasta_file.write(f'>{row["SeqID"]}\n')
            fasta_file.write(f'{row["NuclSeq"]}\n')
    
    print(f'Created {fasta_filename}')

print('FASTA files creation completed.')
