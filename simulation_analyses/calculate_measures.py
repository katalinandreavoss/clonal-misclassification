import sys
import yaml
import pandas as pd

#path = sys.argv[1] 


path="/panfs/qcb-panasas/kavoss/method_comparison/20/0_2/50/1_0/21"
mixcr = path+"clean.fasta.vdjca.clns_IGH.tsv"  # .yaml file from partis partition as input
changeo = path+"vquest_files/combined_db-pass_clone-pass.tsv"
scoper1 = path+"results_db_idClones.tsv"
scoper2 = path+"results_db_hierClones.tsv"
scoper3 = path+"results_db_specClones.tsv"
real_values = path+"family_sizes.txt"