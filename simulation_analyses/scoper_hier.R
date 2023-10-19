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


get_scoper<-function(db,limit,out) {
  results_hierClones <- hierarchicalClones(db, threshold=limit)
  results_db_hierClones <- as.data.frame(results_hierClones)
  l=gsub("\\.","_",as.character(limit))
  path=paste0(out,"results_hierClones_",l,".tsv")
  write.table(results_db_hierClones, file=path, sep='\t',row.names=FALSE,quote=FALSE)
}
limits<- c(0.05,0.1,0.2,0.25,0.3,0.4,0.5)

for (x in limits) {
  get_scoper(test_db,x,opt$out)
}