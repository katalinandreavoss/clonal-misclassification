import sys
import pandas as pd
import statistics
from Bio import SeqIO
import numpy as np

path = sys.argv[1]
#path = "/Users/kavoss/Documents/Research/test_data/44/"
params = path.split("/")
sim = params[-2]
balance = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]

all=path+"all.fasta"
naive=path+"naive.fasta"


def get_id_seq(filepath):
    with open(filepath) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seq = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            seq.append(str(seq_record.seq[:]))
    return identifiers,seq

def hamming_distance(string1, string2):
    # Start with a distance of zero, and count up
    distance = 0
    # Loop over the indices of the string
    L = len(string1)
    for i in range(L):
        # Add 1 to the distance if these two characters are not equal
        if string1[i] != string2[i]:
            distance += 1
    # Return the final count of differences
    return distance

real_id,real_seq = get_id_seq(naive)
real_data = dict(zip(real_id,real_seq))

with open(all) as fasta_file:  # Will close handle cleanly
    identifiers = []
    seq = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        seq.append(str(seq_record.seq[:]))

df = pd.DataFrame({"id": identifiers,"seq":seq})
df['id'] = df['id'].str.split('_c').str[0]
df['real_seq'] = df['id'].map(real_data)
df['hamming_dist'] = df.apply(lambda x: hamming_distance(x.seq, x.real_seq), axis=1)
mean = np.mean(df['hamming_dist'])

result = pd.DataFrame({"SHM": [SHM],"clones":[clones],"leaves": [leaves], "mean_mutations":[mean]})

result.to_csv(path+"mutations.tsv", sep='\t',index=False)