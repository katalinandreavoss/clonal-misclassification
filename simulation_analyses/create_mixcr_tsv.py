import sys
from Bio import SeqIO
import glob


path = sys.argv[1]
#path = "/Users/kavoss/Documents/Research/test_data/44/"
params = path.split("/")
sim = params[-2]
balance = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]

mixcr_path=path+"mixcr_fastas/"
mixcr = open(path+"mixcr.tsv", 'w')
mixcr.write("sequence_id\tclone_id\n")

def get_id_seq(filepath):
    with open(filepath) as fasta_file:  # Will close handle cleanly
        identifiers = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
    return identifiers

files = glob.glob(mixcr_path+"family*_aligned.fasta")
print(files)
i=1
for file in files:
    ids = get_id_seq(file)
    for id in ids:
        mixcr.write(id+"\t"+str(i)+"\n")
    i = i+1


mixcr.flush()
mixcr.close()
