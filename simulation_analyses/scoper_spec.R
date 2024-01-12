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
test_db<-db[ , !names(db) %in% 
               c('clone_id', 'vj_group', 'vjl_group', 'junction_l', 'cdr3_col', 'clone_temp', 'cell_id_temp')]


results_specClones <- spectralClones(test_db, method="vj", germline="germline_alignment_d_mask")
results_db_specClones <- as.data.frame(results_specClones)
write.table(results_db_specClones, paste0(opt$out,"results_specClones_vj.tsv"), sep='\t',row.names=FALSE,quote=FALSE)
