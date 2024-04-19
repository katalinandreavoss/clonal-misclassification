import sys
import pandas as pd
import statistics
from Bio import SeqIO
from scipy.spatial import distance
from itertools import zip_longest
import enchant 
path = sys.argv[1]
params = path.split("/")
sim = params[-2]
balance = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]

naive=path+"naive.fasta"
correct=path+"ancestral_sequences/root_naive.txt"
correct_parsim=path+"ancestral_sequences/ancestral_seqs_parsim.tsv"
correct_ml=path+"ancestral_sequences/ancestral_seqs_ml.tsv"
correct_midpoint= path+"ancestral_sequences/root_naive_rerooted.txt"
changeo=path+"changeo/root_naive.txt"
changeo_midpoint=path+"changeo/root_naive_rerooted.txt"
mixcr=path+"mixcr_fastas/root_naive.txt"
mixcr_midpoint=path+"mixcr_fastas/root_naive_rerooted.txt"
scoper_hier=path+"scoper_hier/root_naive.txt"
scoper_hier_midpoint=path+"scoper_hier/root_naive_rerooted.txt"
scoper_sp=path+"scoper_sp/root_naive.txt"
scoper_sp_midpoint=path+"scoper_sp/root_naive_rerooted.txt"
mptp=path+"mptp/root_naive.txt"
mptp_midpoint=path+"mptp/root_naive_rerooted.txt"

def get_id_seq(filepath):
    with open(filepath) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seq = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            seq.append(str(seq_record.seq[:]))
    return identifiers,seq

real_id,real_seq = get_id_seq(naive)
real_data = dict(zip(real_id,real_seq))



def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip_longest(s1, s2))

def same_length(s1,s2):
    if len(s1)==len(s2):
        return True
    else:
        return False


def get_df(filepath,tool,midpoint):
    df = pd.read_csv(filepath, sep='\t',names=["family", "seq"])
    df['seq'] = df['seq'].astype('str')
    df['tool'] = tool
    df['family'] = df['family'].str.split('_c').str[0]
    df['real_seq'] = df['family'].map(real_data)
    df['real_seq'] = df['real_seq'].astype('str')
    df['hamming_dist'] = df.apply(lambda x: hamming_distance(x.seq, x.real_seq), axis=1)
    df['lev_dist'] = df.apply(lambda x: enchant.utils.levenshtein(x.seq, x.real_seq), axis=1)
    df['len'] = df['seq'].apply(lambda x: len(x))
    df['seq_similarity'] = (1-(df['hamming_dist']/df['len']))*100
    df['seq_similarity_lev'] = (1-(df['lev_dist']/df['len']))*100
    df['same_length'] = df.apply(lambda x: same_length(x.seq, x.real_seq), axis=1)
    df['midpoint_root'] = midpoint
    return df

correct_df = get_df(correct,"correct",False)
correct_parsim_df = get_df(correct_parsim,"correct_parsim",False)
correct_ml_df = get_df(correct_ml,"correct_ml",False)
correct_midpoint_df = get_df(correct_midpoint,"correct",True)
mixcr_df = get_df(mixcr,"mixcr",False)
mixcr_midpoint_df = get_df(mixcr_midpoint,"mixcr",True)
changeo_df = get_df(changeo,"changeo",False)
changeo_midpoint_df = get_df(changeo_midpoint,"changeo",True)
scoper_hier_df = get_df(scoper_hier,"scoper_hier",False)
scoper_hier_midpoint_df = get_df(scoper_hier_midpoint,"scoper_hier",True)
scoper_sp_df = get_df(scoper_sp,"scoper_sp",False)
scoper_sp_midpoint_df = get_df(scoper_sp_midpoint,"scoper_sp",True)
mptp_df = get_df(mptp,"mptp",False)
mptp_midpoint_df = get_df(mptp_midpoint,"mptp",True)


frames = [correct_df, mixcr_df,changeo_df,scoper_hier_df,scoper_sp_df,correct_parsim_df,correct_ml_df,correct_midpoint_df,mixcr_midpoint_df,changeo_midpoint_df,scoper_hier_midpoint_df,scoper_sp_midpoint_df,mptp_df,mptp_midpoint_df]

result = pd.concat(frames)
result = result[result.family != "*"]

result['clones'] = clones
result['SHM'] = SHM
result['leaves'] = leaves


result.to_csv(path+"seq_similarity.tsv", sep='\t',index=False)