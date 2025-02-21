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

# Function to compute overall metrics
def compute_overall_metrics(df):
    clone_family_combos = df.groupby(['clone_id', 'family']).size().reset_index(name='count')
    
    
    # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), True Negatives (TN)
    tp = clone_family_combos['count'].sum()
    fp = (df['clone_id'].value_counts() - df.groupby('clone_id')['family'].count()).sum()
    fn = (df['family'].value_counts() - df.groupby('family')['clone_id'].count()).sum()
    tn = (df.shape[0] - tp - fp - fn)

    # Sensitivity, Precision, F1 Score
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0
    rand_index = (tp + tn) / (tp + fp + fn + tn) if (tp + fp + fn + tn) > 0 else 0

    # Adjusted Rand Index
    ari = adjusted_rand_score(df['clone_id'], df['family'])

    return {
        "sensitivity": sensitivity,
        "precision": precision,
        "f1": f1,
        "rand_index": rand_index,
        "ari": ari
    }

# Compute metrics for all data
metrics_all = compute_overall_metrics(df)
print(f"Metrics with all data: {metrics_all}")

# Remove rows where either clone_id or family is a singleton
clone_grouped = df.groupby(['clone_id'], as_index=False)['sequence_id'].count()
family_grouped = df.groupby(['family'], as_index=False)['sequence_id'].count()

non_singleton_clones = clone_grouped[clone_grouped['sequence_id'] > 1]['clone_id']
non_singleton_families = family_grouped[family_grouped['sequence_id'] > 1]['family']

df_without_singletons = df.loc[df['clone_id'].isin(non_singleton_clones) & df['family'].isin(non_singleton_families)].copy()

# Compute metrics for data without singletons
metrics_no_singletons = compute_overall_metrics(df_without_singletons)
print(f"Metrics without singletons: {metrics_no_singletons}")

# Save results or continue with further processing
print("Done.")
