import sys
import pandas as pd
import statistics
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

changeo = pd.read_csv(changeo_path,sep='\t')
changeo_num_fam = max(changeo["clone_id"])
changeo_grouped = changeo.groupby(['clone_id'])['sequence_id'].count()
changeo_med_size,changeo_var_size = get_med_var(changeo_grouped)
changeo_singletons = changeo.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) == 1).count()
complete.write("changeo\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(changeo_num_fam)+"\t"+str(changeo_med_size)+"\t"+str(changeo_var_size)+"\t"+str(changeo_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

#without singletons
changeo_num_fam = changeo.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) > 1).count()
changeo_grouped = changeo_grouped[changeo.groupby(['clone_id'])['sequence_id'].apply(lambda x: len(x) > 1)]
changeo_med_size,changeo_var_size = get_med_var(changeo_grouped)
without_singletons.write("changeo\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(changeo_num_fam)+"\t"+str(changeo_med_size)+"\t"+str(changeo_var_size)+"\t"+str(changeo_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


scoper1 = pd.read_csv(scoper1_path,sep='\t')
scoper1_num_fam = max(scoper1["clone_id"])
scoper1_grouped = scoper1.groupby(['clone_id'])['sequence_id'].count()
scoper1_med_size,scoper1_var_size = get_med_var(scoper1_grouped)
scoper1_singletons = scoper1.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) == 1).count()
complete.write("scoper_id\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper1_num_fam)+"\t"+str(scoper1_med_size)+"\t"+str(scoper1_var_size)+"\t"+str(scoper1_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

#without singletons
scoper1_num_fam = scoper1.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) > 1).count()
scoper1_grouped = scoper1_grouped[scoper1.groupby(['clone_id'])['sequence_id'].apply(lambda x: len(x) > 1)]
scoper1_med_size,scoper1_var_size = get_med_var(scoper1_grouped)
without_singletons.write("scoper_id\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper1_num_fam)+"\t"+str(scoper1_med_size)+"\t"+str(scoper1_var_size)+"\t"+str(scoper1_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


scoper2 = pd.read_csv(scoper2_path,sep='\t')
scoper2_num_fam = max(scoper2["clone_id"])
scoper2_grouped = scoper2.groupby(['clone_id'])['sequence_id'].count()
scoper2_med_size, scoper2_var_size = get_med_var(scoper2_grouped)
scoper2_singletons = scoper2.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) == 1).count()
complete.write("scoper_hier\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper2_num_fam)+"\t"+str(scoper2_med_size)+"\t"+str(scoper2_var_size)+"\t"+str(scoper2_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

#without singletons
scoper2_num_fam = scoper2.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) > 1).count()
scoper2_grouped = scoper2_grouped[scoper2.groupby(['clone_id'])['sequence_id'].apply(lambda x: len(x) > 1)]
scoper2_med_size, scoper2_var_size = get_med_var(scoper2_grouped)
without_singletons.write("scoper_hier\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper2_num_fam)+"\t"+str(scoper2_med_size)+"\t"+str(scoper2_var_size)+"\t"+str(scoper2_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


scoper3 = pd.read_csv(scoper3_path,sep='\t')
scoper3_num_fam = max(scoper3["clone_id"])
scoper3_grouped = scoper3.groupby(['clone_id'])['sequence_id'].count()
scoper3_med_size, scoper3_var_size = get_med_var(scoper3_grouped)
scoper3_singletons = scoper3.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) == 1).count()
complete.write("scoper_spec\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper3_num_fam)+"\t"+str(scoper3_med_size)+"\t"+str(scoper3_var_size)+"\t"+str(scoper3_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

#without singletons
scoper3_num_fam = scoper3.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) > 1).count()
scoper3_grouped = scoper3_grouped[scoper3.groupby(['clone_id'])['sequence_id'].apply(lambda x: len(x) > 1)]
scoper3_med_size, scoper3_var_size = get_med_var(scoper3_grouped)
without_singletons.write("scoper_spec\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper3_num_fam)+"\t"+str(scoper3_med_size)+"\t"+str(scoper3_var_size)+"\t"+str(scoper3_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

scoper_hier = ["0_05","0_1","0_2","0_25","0_3","0_4","0_5"]
for scop in scoper_hier:
    print(scop)
    scop_path = path+"results_hierClones_"+scop+".tsv"

    scoper = pd.read_csv(scop_path,sep='\t')
    scoper_num_fam = max(scoper["clone_id"])
    scoper_grouped = scoper.groupby(['clone_id'])['sequence_id'].count()
    scoper_med_size, scoper_var_size = get_med_var(scoper_grouped)
    scoper_singletons = scoper.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) == 1).count()
    complete.write("scoper_"+scop+"\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper_num_fam)+"\t"+str(scoper_med_size)+"\t"+str(scoper_var_size)+"\t"+str(scoper_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")

    #without singletons
    scoper_num_fam = scoper.groupby(['clone_id'])['sequence_id'].filter(lambda x: len(x) > 1).count()
    scoper_grouped = scoper_grouped[scoper.groupby(['clone_id'])['sequence_id'].apply(lambda x: len(x) > 1)]
    scoper_med_size, scoper_var_size = get_med_var(scoper_grouped)
    without_singletons.write("scoper_"+scop+"\t"+clones+"\t"+SHM+"\t"+leaves+"\t"+balance+"\t"+str(scoper_num_fam)+"\t"+str(scoper_med_size)+"\t"+str(scoper_var_size)+"\t"+str(scoper_singletons)+"\t"+str(real_values_num_fam)+"\t"+str(real_values_med_size)+"\t"+str(real_values_var_size)+"\n")


complete.flush()
complete.close()
without_singletons.flush()
without_singletons.close()
