import sys
import pandas as pd
from Bio import SeqIO
import numpy as np
from itertools import zip_longest
#import enchant
from concurrent.futures import ThreadPoolExecutor

path = sys.argv[1]
#path = "/Users/kavoss/Documents/Research/14/"
params = path.split("/")
sim = params[-2]
junction_length = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]

naive= path+"naive.fasta"
fasta=path+"clean.fasta"

def get_id_seq(filepath):
    with open(filepath) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seq = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            seq.append(str(seq_record.seq[:]))
    return identifiers,seq


def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip_longest(s1, s2))


def calculate_sim(naive_df, row):
    curr_seq=row['seq']
    curr_fam=row['family_id']
    naive_seq=naive_df.loc[naive_df['id'] == curr_fam, 'seq'].item()
    print(curr_seq)
    print(naive_seq)
    ham = hamming_distance(naive_seq,curr_seq)
    return  ham/len(curr_seq)*100


real_id,real_seq = get_id_seq(fasta)

real_data = pd.DataFrame(
    {'id': real_id,
     'seq': real_seq
    })
real_data[['family_id','clone_id']] = real_data['id'].str.split('_c',expand=True)

naive_id,naive_seq = get_id_seq(naive)

naive_data = pd.DataFrame(
    {'id': naive_id,
     'seq': naive_seq
    })


real_data['distance'] = real_data.apply(lambda row: calculate_sim(naive_data, row), axis=1)


print(real_data)
real_data['SHM'] = SHM
real_data['leaves'] = leaves
real_data = real_data.drop("clone_id", axis='columns')
real_data = real_data.drop("seq", axis='columns')
real_data = real_data.drop("id", axis='columns')
real_data = real_data.drop("family_id", axis='columns')

real_data.to_csv(path+"distance.tsv", sep='\t',index=False)

