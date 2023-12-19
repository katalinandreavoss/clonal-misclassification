import sys
import pandas as pd
import statistics
from Bio import SeqIO


path = sys.argv[1]
#path = "/Users/kavoss/Documents/Research/simulations/14/0_2/10/6/"
params = path.split("/")
sim = params[-2]
balance = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]


with open(path+"clean.fasta") as fasta_file:  # Will close handle cleanly
    identifiers = []
    seq = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        seq.append(str(seq_record.seq[:]))

fasta = pd.DataFrame({"id": identifiers,"seq":seq})

changeo_path = path+"vquest_files/combined_db-pass_clone-pass.tsv"
scoper2_path = path+"results_hierClones.tsv"
scoper3_path = path+"results_specClones.tsv"
mptp_path = path+"mptp_data_singletons.txt"

all = [changeo_path,scoper2_path,scoper3_path,mptp_path]
file_path=pd.DataFrame({"file": all,"dir":["changeo","scoper_hier","scoper_sp","mptp"]})



for index, row in file_path.iterrows():
    df = pd.read_csv(row["file"],sep='\t')
    df_grouped = df.groupby(['clone_id'], as_index = False)['sequence_id'].count()
    clones=df_grouped[df_grouped["sequence_id"] != 1]["clone_id"]
    for clone in clones:
        ids=df[df["clone_id"]==clone]["sequence_id"]
        file = open(path + row["dir"]+"/" + ids.iloc[1]+".tsv", 'w')
        for id in ids:
            file.write(">"+str(id)+"\n"+str(fasta[fasta["id"]==id]["seq"].values[0]+"\n"))
        file.flush()
        file.close()





