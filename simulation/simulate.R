#!/usr/bin/env Rscript
library("optparse")
library("Biostrings")
library(tidyr)
library(TreeSim)
library(ape)
require(phytools)
require(splits)
require(gtools)
library(phangorn)

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="path to directory", metavar="character"),
  make_option(c("-r", "--rate"), type="character", 
              help="SHM rate", metavar="character"),
  make_option(c("-p", "--path"), type="character", 
              help="path to VDJ directory containing fastas with VDJ genes", metavar="character"),
  make_option(c("-l", "--leaves"), type="integer", 
              help="number of leaves", metavar="character"),
  make_option(c("-j", "--junction"), type="integer", 
              help="length of junction region", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


V_file <- readDNAStringSet(paste0(opt$path,"/IGHV.fasta"))
seq_name = names(V_file)
sequence = paste(V_file)
V_genes <- data.frame(seq_name, sequence)
V_genes$sequence<-gsub("\\.","",as.character(V_genes$sequence))
D_file <- readDNAStringSet(paste0(opt$path,"/IGHD.fasta"))
seq_name = names(D_file)
sequence = paste(D_file)
D_genes <- data.frame(seq_name, sequence)
D_genes$sequence<-gsub("\\.","",D_genes$sequence)
J_file <- readDNAStringSet(paste0(opt$path,"/IGHJ.fasta"))
seq_name = names(J_file)
sequence = paste(J_file)
J_genes <- data.frame(seq_name, sequence)
J_genes$sequence<-gsub("\\.","",J_genes$sequence)


filepath<-paste0(opt$dir,"/clean.fasta")
file <- file(description = filepath, open = "a")

generate_random_dna_sequence <- function(length) {
  bases <- c("a", "t", "c", "g")
  random_sequence <- paste(sample(bases, length, replace = TRUE), collapse = "")
  return(random_sequence)
}

for (x in 1:10) {
  V<-sample(V_genes$sequence,1)
  D<-sample(D_genes$sequence,1)
  J<-sample(J_genes$sequence,1)

  leaves<-opt$leaves
  
  trees<-sim.bd.taxa(n=leaves,numbsim=1,lambda=2.0,mu=0.5)
  tree<-trees[[1]]
  rate<-as.numeric(gsub("_",".",opt$rate))
  junction_rest<-as.numeric(opt$junction)
  j_half<-floor(junction_rest/2)

 
  junction <- generate_random_dna_sequence(junction_rest)
  junction_region <- paste0(substring(junction,1,j_half),tolower(D),substring(junction,j_half+1,junction_rest))
  
  junction_r_char<-unlist(strsplit(junction_region, ""))

  sims <- simSeq(tree, l=length(junction_r_char), type="DNA",
                  rootseq = junction_r_char, rate=rate)
  
  filepath_fam<-paste0(opt$dir,"/family_",x,".fasta")
  file_fam <- file(description = filepath_fam, open = "a")

  for (s in 1:length(sims)){
    curr_junction<-toupper(paste(as.character(sims[s]), collapse = ''))
    sequence<-paste0(V,curr_junction,J)
    cat(paste0(">family_",x,"_clone_",s), file= file, append=TRUE,sep="\n")
    cat(toupper(sequence), file= file, append=TRUE,sep="\n")
    cat(paste0(">family_",x,"_clone_",s), file= file_fam, append=TRUE,sep="\n")
    cat(toupper(sequence), file= file_fam, append=TRUE,sep="\n")
  }
  close(file_fam)
}

close(file)

