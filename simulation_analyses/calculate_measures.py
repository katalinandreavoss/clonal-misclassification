import sys
import pandas as pd
import statistics
import os
path = sys.argv[1]

#path = "/Users/kavoss/Documents/Research/simulations/14/0_2/10/6/"
params = path.split("/")
sim = params[-2]
balance = params[-3]
leaves = params[-4]
SHM = params[-5]
clones = params[-6]
complete = open(path+"analysis_complete.tsv", 'w')
without_singletons = open(path+"analysis_no_singletons.tsv", 'w')

complete.write("tool\tclones\tSHM\tleaves\tbalance\tnum_families\tmedian_size_families\tvar_size_families\tsingletons\treal_num_fam\treal_med_size\treal_var_size\n")
without_singletons.write("tool\tclones\tSHM\tleaves\tbalance\tnum_families\tmedian_size_families\tvar_size_families\tsingletons\treal_num_fam\treal_med_size\treal_var_size\n")

mixcr_path = path+"clean.fasta.vdjca.clns_IGH.tsv"  # .yaml file from partis partition as input
changeo_path = path+"vquest_files/combined_db-pass_clone-pass.tsv"
scoper1_path = path+"results_db_idClones.tsv"
scoper2_path = path+"results_hierClones.tsv"
scoper3_path = path+"results_specClones.tsv"
scoper4_path = path+"results_specClones_vj.tsv"
mptp_path= path+"mptp_data_singletons.txt"
gmyc_path= path+"gmyc.tsv"
real_values_path = path+"family_sizes.txt"


real_values = pd.read_csv(real_values_path, header=None)
real_values_num_fam = len(real_values)
real_values_med_size = statistics.median(real_values[4])
real_values_var_size = statistics.variance(real_values[4])

def mixcr_get_med_var(df):
    if len(df)>1:
        mixcr_med_size = statistics.median(df)
        mixcr_var_size = statistics.variance(df)
    elif len(df)==1:
        mixcr_med_size = statistics.median(df)
        mixcr_var_size = "NaN"
    else:
        mixcr_med_size = "NaN"
        mixcr_var_size = "NaN"
    return mixcr_med_size, mixcr_var_size

mixcr = pd.read_csv(mixcr_path,sep='\t')
mixcr_num_fam = len(mixcr)
mixcr_med_size, mixcr_var_size = mixcr_get_med_var(mixcr["readCount"])
mixcr_singletons = len(mixcr[mixcr["readCount"]==1])
complete.write("mixcr\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(mixcr_num_fam)+"\t"+str(mixcr_med_size)+"\t"+str(mixcr_var_size)+"\t"+str(mixcr_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

#without singletons
mixcr_num_fam = len(mixcr[mixcr["readCount"]!=1])
mixcr_singletons = len(mixcr[mixcr["readCount"]!=1])
mixcr_med_size, mixcr_var_size = mixcr_get_med_var(mixcr[mixcr["readCount"]!=1]["readCount"])
without_singletons.write("mixcr\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(mixcr_num_fam)+"\t"+str(mixcr_med_size)+"\t"+str(mixcr_var_size)+"\t"+str(mixcr_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


def get_med_var(df):
    if len(df)>1:
        med_size = statistics.median(df)
        var_size = statistics.variance(df)
    elif len(df)==1:
        med_size = statistics.median(df)
        var_size = "NaN"
    else:
        med_size = "NaN"
        var_size = "NaN"
    return med_size,var_size

def get_results_all(path):
    df = pd.read_csv(path, sep='\t')
    num_fam = max(df["clone_id"])
    df_grouped = df.groupby(['clone_id'])['sequence_id'].count()
    med_size, var_size = get_med_var(df_grouped)
    singletons = df.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) == 1).count()
    # without singletons
    non_singletons = df.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) > 1)
    num_fam_no_s = len(pd.unique(df.loc[df['sequence_id'].isin(non_singletons)]['clone_id']))

    df_grouped = df_grouped[df.groupby(['clone_id'])['sequence_id'].apply(lambda x: len(x) > 1)]
    med_size_no_s, var_size_no_s = get_med_var(df_grouped)

    return num_fam, med_size, var_size, singletons,num_fam_no_s,med_size_no_s, var_size_no_s

changeo_num_fam,changeo_med_size,changeo_var_size,changeo_singletons,changeo_num_fam_no_s,changeo_med_size_no_s,changeo_var_size_no_s = get_results_all(changeo_path)
complete.write("changeo\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(changeo_num_fam)+"\t"+str(changeo_med_size)+"\t"+str(changeo_var_size)+"\t"+str(changeo_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
without_singletons.write("changeo\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(changeo_num_fam_no_s)+"\t"+str(changeo_med_size_no_s)+"\t"+str(changeo_var_size_no_s)+"\t"+str(changeo_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(scoper1_path)
complete.write("scoper_id\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
without_singletons.write("scoper_id\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(scoper2_path)
complete.write("scoper_hier\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
without_singletons.write("scoper_hier\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(scoper3_path)
complete.write("scoper_spec\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
without_singletons.write("scoper_spec\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

# num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(scoper4_path)
# complete.write("scoper_spec_vj\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
# without_singletons.write("scoper_spec_vj\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(mptp_path)
complete.write("mptp\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
without_singletons.write("mptp\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


# if os.path.exists(gmyc_path):
#     num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(gmyc_path)
#     complete.write("gmyc\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
#     without_singletons.write("gmyc\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


# scoper_hier = ["0_05","0_1","0_2","0_25","0_3","0_4","0_5"]
# for scop in scoper_hier:
#     print(scop)
#     scop_path = path+"results_hierClones_"+scop+".tsv"
#     num_fam,med_size,var_size,singletons,num_fam_no_s,med_size_no_s,var_size_no_s = get_results_all(scop_path)
#     complete.write("scoper_"+scop+"\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam)+"\t"+str(med_size)+"\t"+str(var_size)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")
#     without_singletons.write("scoper_"+scop+"\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(num_fam_no_s)+"\t"+str(med_size_no_s)+"\t"+str(var_size_no_s)+"\t"+str(singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


complete.flush()
complete.close()
without_singletons.flush()
without_singletons.close()
