#!/usr/bin/env Rscript
library("optparse")
library(phangorn)
library(ape)
option_list = list(
  make_option(c("-d", "--dir"), type="character", default="out.txt", 
              help="directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path<-opt$dir

filenames <- list.files(path=path, pattern="*_tree_.raxml.bestTree$", full.names=TRUE,recursive = TRUE)

for(file in filenames){
  tree<-read.tree(file)
  family<-gsub(path,"",gsub("_tree_.raxml.bestTree","",file))
  tryCatch({
    midpoint_tree<-midpoint(tree)
    write.tree(midpoint_tree, paste0(path,family,"_tree_rerooted.nexus"))
  },
  error = function(e) {
    print(file)
  }
  )
}

