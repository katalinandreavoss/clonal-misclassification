import sys
import pandas as pd
import os
from Bio import SeqIO

path = sys.argv[1]
#path = "/Users/kavoss/Documents/Research/test_data/44/"
params = path.split("/")
sim = params[-2]
balance = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]

real_values_path = path + "clean.fasta"
mixcr_path = path + "mixcr.tsv"
changeo_path = path + "vquest_files/combined_db-pass_clone-pass.tsv"
scoper_hier_path = path + "results_hierClones.tsv"
scoper_sp_path = path + "results_specClones.tsv"
mptp_path = path + "mptp_data_singletons.txt"


def read_fasta(file_path):
    sequence_names = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_names.append(record.id)
    return sequence_names


sequence_names = read_fasta(real_values_path)
print(sequence_names)

def add_missing_sequences(df, sequence_names):
    max_clone_id = df['clone_id'].max() if not df.empty else 0

    for name in sequence_names:
        if name not in df['sequence_id'].values:
            max_clone_id += 1
            df = pd.concat([df, pd.DataFrame({'sequence_id': name, 'clone_id': max_clone_id}, index=[0])], ignore_index=True)

    return df


all = [mixcr_path, changeo_path, scoper_hier_path, scoper_sp_path, mptp_path]
file_path = pd.DataFrame({"file": all, "tool": ["mixcr", "changeo", "scoper_hier", "scoper_sp", "mptp"]})
file_path["TP"] = 0.0
file_path["TN"] = 0.0
file_path["FN"] = 0.0
file_path["FP"] = 0.0
file_path["sensitivity"] = 0.0
file_path["precision"] = 0.0
file_path["f1"] = 0.0
file_path['SHM'] = SHM
file_path['clones'] = clones
file_path['leaves'] = leaves
file_path['junction_length'] = balance
file_path['sim'] = sim
file_path['junction_length_scoper'] = 0.0
file_path['solved'] = 0.0
file_path['all'] = len(sequence_names)
file_path['unsolved'] = len(sequence_names)
file_path['ratio'] = 0.0
file_path['rand'] = 0.0

junction_length = 0.0

for i, r in file_path.iterrows():
    if os.path.exists(r["file"]):
        df = pd.read_csv(r["file"], sep='\t')
        solved_sequences = len(df)
        df = add_missing_sequences(df, sequence_names)
        all_sequences = len(df)
        df["family"] = df["sequence_id"].str.split('_c').str[0]
        df_grouped = df.groupby(['clone_id'], as_index=False)['sequence_id'].count()
        clones = df_grouped[df_grouped["sequence_id"] != 1]["clone_id"].values
        df = df.loc[df['clone_id'].isin(clones)]
        df["TP"] = 0.0
        df["TN"] = 0.0
        df["FN"] = 0.0
        df["FP"] = 0.0
        df["sensitivity"] = 0.0
        df["precision"] = 0.0
        df["f1"] = 0.0
        df["rand"] = 0.0
        if "hierClones" in r["file"]:
            junction_length = df["junction_length"].mean()
        for index, row in df.iterrows():
            family = row["family"]
            clone = row["clone_id"]
            TP = df.loc[(df["clone_id"] == clone) & (df["family"] == family)].shape[0]
            TN = df.loc[(df["clone_id"] != clone) & (df["family"] != family)].shape[0]
            FN = df.loc[(df["clone_id"] != clone) & (df["family"] == family)].shape[0]
            FP = df.loc[(df["clone_id"] == clone) & (df["family"] != family)].shape[0]
            df.loc[index, 'TP'] = TP
            df.loc[index, 'TN'] = TN
            df.loc[index, 'FN'] = FN
            df.loc[index, 'FP'] = FP
            df.loc[index, 'sensitivity'] = TP / (TP + FN)
            df.loc[index, 'precision'] = TP / (TP + FP)
            df.loc[index, 'f1'] = 2 * TP / (2 * TP + FP + FN)
            df.loc[index, 'rand'] = (TP + TN)/(TP + FP + FN + TN)
        file_path.loc[i, 'TP'] = df["TP"].mean()
        file_path.loc[i, 'TN'] = df["TN"].mean()
        file_path.loc[i, 'FN'] = df["FN"].mean()
        file_path.loc[i, 'FP'] = df["FP"].mean()
        file_path.loc[i, 'precision'] = df["precision"].mean()
        file_path.loc[i, 'sensitivity'] = df["sensitivity"].mean()
        file_path.loc[i, 'f1'] = df["f1"].mean()
        file_path.loc[i, 'rand'] = df["rand"].mean()
        file_path.loc[i, 'solved'] = solved_sequences
        file_path.loc[i, 'all'] = all_sequences
        file_path.loc[i, 'unsolved'] = all_sequences-solved_sequences
        file_path.loc[i, 'ratio'] = solved_sequences/all_sequences

file_path["junction_length_scoper"] = junction_length
file_path = file_path.drop(columns=['file'])
file_path.to_csv(path + "f1_all.tsv", sep='\t', index=False)
