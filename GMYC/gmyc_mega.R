#!/usr/bin/env Rscript
library("optparse")
require(ape)
require(phytools)
require(splits)
require(gtools)
library(phangorn)

option_list = list(
  make_option(c("-d", "--dir"), type="character", default="out.txt", 
              help="directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path<-opt$dir

tree_path<-paste0(path,"tree_files/mega_tree_.raxml.bestTree")
tree<-read.tree(tree_path)

#consensus_ultra=chronos(tree, lambda=0)
consensus_ultra=force.ultrametric(tree, method="extend")

midpoint_tree<-midpoint(consensus_ultra)

single<-gmyc(midpoint_tree,method = "single",quiet = TRUE)
#plot.gmyc(single)

#mulit run really long and is only experimental
#multi<-gmyc(midpoint_tree,method = "m",quiet = TRUE)

specs<-spec.list(single)
colnames(specs)<-c("clone_id","sequence_id")
write.table(specs, file = paste0(path,"gmyc_force.tsv"), row.names=FALSE, sep="\t",quote = FALSE)
