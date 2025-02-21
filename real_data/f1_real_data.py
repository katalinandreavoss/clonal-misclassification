import sys
import pandas as pd
import os

#path = sys.argv[1]

path = "/project/mpennell_978/kavoss/real_data_cf/14007/3/"
params = path.split("/")
timepoint = params[-2]

og_path = path+"mapping.csv"
mptp_path = path+"mptp_data_singletons.txt"


mptp = pd.read_csv(mptp_path,sep='\t')
og = pd.read_csv(og_path)

og.rename(columns={'clone_id': 'family'}, inplace=True)

# Merge the two DataFrames on 'sequence_id'
df = pd.merge(mptp, og, on='sequence_id', how='inner')

df_grouped = df.groupby(['clone_id'], as_index=False)['sequence_id'].count()
clones = df_grouped[df_grouped["sequence_id"] != 1]["clone_id"].values
df = df.loc[df['clone_id'].isin(clones)]
df["TP"]= 0.0
df["TN"] = 0.0
df["FN"] = 0.0
df["FP"] = 0.0
df["sensitivity"] = 0.0
df["precision"] = 0.0
df["f1"] = 0.0
df["rand"] = 0.0
for index, row in df.iterrows():
    family= row["family"]
    clone = row["clone_id"]
    TP = df.loc[(df["clone_id"] == clone) & (df["family"] == family)].shape[0]
    TN = df.loc[(df["clone_id"] != clone) & (df["family"] != family)].shape[0]
    FN = df.loc[(df["clone_id"] != clone) & (df["family"] == family)].shape[0]
    FP = df.loc[(df["clone_id"] == clone) & (df["family"] != family)].shape[0]
    df.loc[index, 'TP'] = TP
    df.loc[index, 'TN'] = TN
    df.loc[index, 'FN'] = FN
    df.loc[index, 'FP'] = FP
    df.loc[index, 'sensitivity'] = TP/(TP+FN)
    df.loc[index, 'precision'] = TP / (TP + FP)
    df.loc[index, 'f1'] = 2*TP / (2*TP + FP + FN)
    df.loc[index, 'rand'] = TP+TN / (TP + FP + FN+TN)


df.to_csv(path+"sensitivity_precision_slow.tsv", sep='\t',index=False)

