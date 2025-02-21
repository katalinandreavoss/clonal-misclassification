import pandas as pd
from sklearn.metrics import adjusted_rand_score

# Define paths
path = "/project/mpennell_978/kavoss/real_data_cf/14007/3/"
params = path.split("/")
timepoint = params[-2]

og_path = path + "mapping.csv"
mptp_path = path + "mptp_data_singletons.txt"

# Load data
mptp = pd.read_csv(mptp_path, sep='\t')
og = pd.read_csv(og_path)

# Rename 'clone_id' column in og to 'family'
og.rename(columns={'clone_id': 'family'}, inplace=True)

# Merge the two DataFrames on 'sequence_id'
df = pd.merge(mptp, og, on='sequence_id', how='inner')

# Function to compute TP, TN, FP, FN and metrics
def compute_metrics(df):
    df["TP"] = 0.0
    df["TN"] = 0.0
    df["FN"] = 0.0
    df["FP"] = 0.0

    clone_family_combos = df.groupby(['clone_id', 'family']).size().reset_index(name='count')
    clone_sizes = df['clone_id'].value_counts()
    family_sizes = df['family'].value_counts()

    total_rows = len(df)
    progress_interval = max(1, total_rows // 10)

    for idx, row in df.iterrows():
        clone = row['clone_id']
        family = row['family']

        TP = clone_family_combos.loc[(clone_family_combos['clone_id'] == clone) & 
                                     (clone_family_combos['family'] == family), 'count'].values[0]
        FN = family_sizes[family] - TP
        FP = clone_sizes[clone] - TP
        TN = df.shape[0] - TP - FN - FP

        df.loc[idx, 'TP'] = TP
        df.loc[idx, 'TN'] = TN
        df.loc[idx, 'FN'] = FN
        df.loc[idx, 'FP'] = FP

        if (idx + 1) % progress_interval == 0 or idx + 1 == total_rows:
            progress = min(100, ((idx + 1) / total_rows) * 100)
            print(f"Progress: {progress:.2f}%")

    df["sensitivity"] = df["TP"] / (df["TP"] + df["FN"])
    df["precision"] = df["TP"] / (df["TP"] + df["FP"])
    df["f1"] = 2 * df["TP"] / (2 * df["TP"] + df["FP"] + df["FN"])
    df["rand"] = (df["TP"] + df["TN"]) / (df["TP"] + df["FP"] + df["FN"] + df["TN"])

    return df

# Process with singletons
print("Processing with singletons:")
df_with_singletons = df.copy()
#df_with_singletons = compute_metrics(df_with_singletons)

# Remove rows where either clone_id or family is a singleton
clone_grouped = df.groupby(['clone_id'], as_index=False)['sequence_id'].count()
family_grouped = df.groupby(['family'], as_index=False)['sequence_id'].count()

non_singleton_clones = clone_grouped[clone_grouped['sequence_id'] > 1]['clone_id']
non_singleton_families = family_grouped[family_grouped['sequence_id'] > 1]['family']

print("Processing without singletons:")
df_without_singletons = df.loc[df['clone_id'].isin(non_singleton_clones) & df['family'].isin(non_singleton_families)].copy()
df_without_singletons = compute_metrics(df_without_singletons)

#df_with_singletons.to_csv(path+"sensitivity_precision.tsv", sep='\t',index=False)
df_without_singletons.to_csv(path+"sensitivity_precision_no_singletons.tsv", sep='\t',index=False)
print("Done.")