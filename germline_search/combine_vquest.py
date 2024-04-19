import sys
from glob import glob
import pandas as pd
import numpy as np
path = sys.argv[1]  #vquest directory
subfolders=glob(path+"/*/", recursive = True)

first1 = pd.read_table(path+"/001/1_Summary.txt", delimiter='\t',index_col=0)
first2 = pd.read_table(path+"/001/2_IMGT-gapped-nt-sequences.txt", delimiter='\t',index_col=0)
first3 = pd.read_table(path+"/001/3_Nt-sequences.txt", delimiter='\t',index_col=0)
first6 = pd.read_table(path+"/001/6_Junction.txt", delimiter='\t',index_col=0)

first1["V-REGION identity % (with ins/del events)"] = ""
first1["V-REGION identity nt (with ins/del events)"] = ""
first1["V-REGION insertions"] = ""
first1["V-REGION deletions"] = ""
first6 = first6.fillna(0)
first6 = first6.apply(lambda x: x.astype(int) if x.dtype == 'float' else x)


for sub in subfolders:
   if "combined" not in sub and "/001/" not in sub:
        file1 =pd.read_table(sub + "1_Summary.txt", delimiter='\t', index_col=0)
        file1["V-REGION identity % (with ins/del events)"] = ""
        file1["V-REGION identity nt (with ins/del events)"] = ""
        file1["V-REGION insertions"] = ""
        file1["V-REGION deletions"] = ""
        frames = [first1, file1]
        first1 = pd.concat(frames)
        first1 = first1.loc[:, ~first1.columns.str.contains('^Unnamed')]
        first1 = first1.reset_index(drop=True)
        first1 = first1.rename_axis('Sequence Number')
        first1.index += 1

        file2 = pd.read_table(sub + "2_IMGT-gapped-nt-sequences.txt", delimiter='\t', index_col=0)
        frames = [first2, file2]
        first2 = pd.concat(frames)
        first2 = first2.loc[:, ~first2.columns.str.contains('^Unnamed')]
        first2 = first2.reset_index(drop=True)
        first2 = first2.rename_axis('Sequence Number')
        first2.index += 1

        file3 = pd.read_table(sub + "3_Nt-sequences.txt", delimiter='\t', index_col=0)
        frames = [first3, file3]
        first3 = pd.concat(frames)
        first3 = first3.loc[:, ~first3.columns.str.contains('^Unnamed')]
        first3 = first3.reset_index(drop=True)
        first3 = first3.rename_axis('Sequence Number')
        first3.index += 1

        file6 = pd.read_table(sub + "6_Junction.txt", delimiter='\t', index_col=0)
        frames = [first6, file6]
        first6 = pd.concat(frames)
        first6 = first6.loc[:, ~first6.columns.str.contains('^Unnamed')]
        first6 = first6.reset_index(drop=True)
        first6 = first6.rename_axis('Sequence Number')
        first6.index += 1
        first6 = first6.fillna(0)
        first6 = first6.apply(lambda x: x.astype(int) if x.dtype == 'float' else x)

first3 = first3.fillna(0)
first3['V-REGION start'] = first3['V-REGION start'].apply(np.int64)

first1.to_csv(path+"/combined/1_Summary.txt", sep="\t")
first2.to_csv(path+"/combined/2_IMGT-gapped-nt-sequences.txt", sep="\t")
first3.to_csv(path+"/combined/3_Nt-sequences.txt", sep="\t")
first6.to_csv(path+"/combined/6_Junction.txt", sep="\t")
