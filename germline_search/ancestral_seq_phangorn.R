#!/usr/bin/env Rscript
library("optparse")
library(phangorn)
library(ape)
library("Biostrings")
library(tidyverse)
option_list = list(
  make_option(c("-d", "--dir"), type="character", default="out.txt", 
              help="directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path<-opt$dir

filenames <- list.files(path=paste0(path,"tree_files"), pattern="*_tree_.raxml.bestTree$", full.names=TRUE,recursive = TRUE)

for(file in filenames){
  tree<-read.tree(file)
  family<-gsub(paste0(path,"tree_files/"),"",gsub("_tree_.raxml.bestTree","",file))
  data<-read.phyDat(paste0(path,family,".fasta"), format = "fasta", type = "DNA")
  anc.mpr <- ancestral.pars(tree, data, "MPR", return = "phyDat")
  write.phyDat(anc.mpr, file = paste0(path,family,"_ancestral_parsimony.fasta"), format = "fasta")
  fit <- pml(tree,data)
  fit <- optim.pml(fit, model="F81", control = pml.control(trace=0))
  anc.ml <- ancestral.pml(fit, "ml", return = "phyDat")
  write.phyDat(anc.ml, file = paste0(path,family,"_ancestral_ml.fasta"), format = "fasta")
}

parsim <- list.files(path=path, pattern="*_ancestral_parsimony.fasta$", full.names=TRUE,recursive = TRUE)
ml <- list.files(path=path, pattern="*_ancestral_ml.fasta$", full.names=TRUE,recursive = TRUE)

parsim_df = data.frame(family=character(0),sequence=character(0))
ml_df = data.frame(family=character(0),sequence=character(0))


for(file in parsim){
  fastaFile <- readDNAStringSet(file)
  family<-gsub(paste0(path,"/"),"",gsub("_ancestral_parsimony.fasta","",file))
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  nodes<-df[! grepl("family", df$seq_name), ]
  parsim_df<-parsim_df %>% add_row(family = family, sequence = nodes[1,]$sequence)
}

for(file in ml){
  fastaFile <- readDNAStringSet(file)
  family<-gsub(paste0(path,"/"),"",gsub("_ancestral_ml.fasta","",file))
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  nodes<-df[! grepl("family", df$seq_name), ]
  ml_df<-ml_df %>% add_row(family = family, sequence = nodes[1,]$sequence)
}

write.table(parsim_df, file=paste0(path,"ancestral_sequences/ancestral_seqs_parsim.tsv"), quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)
write.table(ml_df, file=paste0(path,"ancestral_sequences/ancestral_seqs_ml.tsv"), quote=FALSE, sep='\t', row.names = FALSE,col.names = FALSE)




