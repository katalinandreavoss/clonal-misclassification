import sys
import pandas as pd
from Bio import SeqIO
import numpy as np
from itertools import zip_longest
#import enchant
from concurrent.futures import ThreadPoolExecutor

path = sys.argv[1]
#path = "/Users/kavoss/Documents/Research/"
params = path.split("/")
sim = params[-2]
junction_length = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]

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


def calculate_min_max(df, curr_id, curr_seq, curr_fam):
    # Filter rows based on family_id

    same_family_rows = df[(df['family_id'] == curr_fam) & (df['id'] != curr_id)].copy()
    diff_family_rows = df[df['family_id'] != curr_fam].copy()


    # Calculate Hamming distance and Levenshtein distance for same_family_rows
    same_family_rows['ham_dist'] = same_family_rows['seq'].apply(lambda x: hamming_distance(curr_seq, x))
    same_family_rows['ham_dist_norm'] = same_family_rows['ham_dist'] / len(curr_seq)
    same_family_rows['ham_dist_norm_j'] = same_family_rows['ham_dist'] / int(junction_length)

    #same_family_rows['lev_dist'] = same_family_rows['seq'].apply(lambda x: enchant.utils.levenshtein(curr_seq, x) / len(curr_seq))

    diff_family_rows['ham_dist'] = diff_family_rows['seq'].apply(lambda x: hamming_distance(curr_seq, x))
    diff_family_rows['ham_dist_norm'] = diff_family_rows['ham_dist'] / len(curr_seq)
    diff_family_rows['ham_dist_norm_j'] = diff_family_rows['ham_dist'] / int(junction_length)

    #diff_family_rows['lev_dist'] = diff_family_rows['seq'].apply(lambda x: enchant.utils.levenshtein(curr_seq, x) / len(curr_seq))


    # Calculate max_within, min_between, max_within_lev, min_between_lev
    max_within_j = same_family_rows['ham_dist_norm_j'].max()
    min_between_j = diff_family_rows['ham_dist_norm_j'].min()
    max_within = same_family_rows['ham_dist_norm'].max()
    min_between = diff_family_rows['ham_dist_norm'].min()
    #max_within_lev = same_family_rows['lev_dist'].max()
    #min_between_lev = diff_family_rows['lev_dist'].min()

    return max_within, min_between,max_within_j,min_between_j#, max_within_lev, min_between_lev


real_id,real_seq = get_id_seq(fasta)

real_data = pd.DataFrame(
    {'id': real_id,
     'seq': real_seq
    })


real_data[['family_id','clone_id']] = real_data['id'].str.split('_c',expand=True)

def process_row(row):
    return calculate_min_max(real_data, row.id, row.seq, row.family_id)


with ThreadPoolExecutor(max_workers=100) as executor:
    results = list(executor.map(process_row, real_data.itertuples(index=False)))





result_all = pd.DataFrame(results, columns=['max_within', 'min_between','max_within_j','min_between_j'])
result_df1 = result_all[['max_within', 'min_between']]
reshaped_df1 = pd.melt(result_df1, var_name='category', value_name='ham_dist')

result_df2 = result_all[['max_within_j','min_between_j']] #'max_within_lev', 'min_between_lev'
reshaped_df2 = pd.melt(result_df2, var_name='category', value_name='ham_dist_j')

reshaped_df = pd.concat([reshaped_df1, reshaped_df2], axis=1)

reshaped_df[["SHM"]] = SHM
reshaped_df[["leaves"]] = leaves
reshaped_df[["junction_length"]] = junction_length
reshaped_df.to_csv(path+"ham_distance.tsv", sep='\t',index=False)

