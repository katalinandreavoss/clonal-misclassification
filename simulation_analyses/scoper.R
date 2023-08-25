#!/usr/bin/env Rscript
library("optparse")
library("scoper")
library(dplyr)

option_list = list(
  make_option(c("-d", "--db"), type="character", default=NULL, 
              help="path to changeo database", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


db<-read.table(opt$db, header=TRUE, fill=TRUE,sep="\t")
db$germline_alignment_d_mask<-NA
test_db<-db[ , !names(db) %in% 
               c('clone_id', 'vj_group', 'vjl_group', 'junction_l', 'cdr3_col', 'clone_temp', 'cell_id_temp')]


results_idClones <- identicalClones(test_db, method="nt")
results_db_idClones <- as.data.frame(results_idClones)
write.table(results_db_idClones, paste0(opt$out,"results_db_idClones.tsv"), sep='\t',row.names=FALSE,quote=FALSE)

results_hierClones <- hierarchicalClones(test_db, threshold=0.15)
results_db_hierClones <- as.data.frame(results_hierClones)
write.table(results_db_hierClones, paste0(opt$out,"results_hierClones.tsv"), sep='\t',row.names=FALSE,quote=FALSE)

results_specClones <- spectralClones(test_db, method="novj")
results_db_specClones <- as.data.frame(results_specClones)
write.table(results_db_specClones, paste0(opt$out,"results_specClones.tsv"), sep='\t',row.names=FALSE,quote=FALSE)
