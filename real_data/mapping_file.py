import pandas as pd
import os

# Function to process SeqID by cutting off anything after ":"
def process_seqid(seqid):
    return seqid.split(":")[0]

# Read the CSV file
input_csv = "/project/mpennell_978/kavoss/real_data_cf/14007/14007.csv"  # Replace with your file path
data = pd.read_csv(input_csv)


output_dir = "/project/mpennell_978/kavoss/real_data_cf/14007/"

# Group the data by TimePoint and create a mapping file for each TimePoint
for timepoint, group in data.groupby('TimePoint'):
    # Process SeqID column
    group['SeqID'] = group['SeqID'].apply(process_seqid)
    
    # Select the required columns
    mapping = group[['SeqID', 'TimePointCloneID']].rename(columns={'SeqID': 'sequence_id', 'TimePointCloneID': 'clone_id'})
    
    # Save the mapping file
    output_file = os.path.join(output_dir, f"mapping_{timepoint}.csv")
    mapping.to_csv(output_file, index=False)

    print(f"Mapping file for TimePoint {timepoint} saved as {output_file}")
