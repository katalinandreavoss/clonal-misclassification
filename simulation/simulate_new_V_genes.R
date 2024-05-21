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
              help="amount of mutations introduced", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt$dir<-"/scratch1/kavoss/sims_fake_V/10/0_01/100/20/8/"
#opt$rate<-"0_01"
#opt$leaves<-100
#opt$path<-"/home1/kavoss/VDJ_data/"
#opt$junction<-"20"
print(opt$junction)

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

mutations<-as.numeric(opt$junction)


mutate_sequence <- function(sequence, mut_avg) {
  nucleotides <- c("A", "C", "G", "T")
  
  # Generate random positions for mutations, insertions, and deletions
  sample_number <- floor(rnorm(1, mean = mut_avg, sd = 5))
  #print(sample_number)
  mutation_positions <- sample.int(nchar(sequence), size = sample_number, replace = FALSE)
  insertion_positions <- sort(sample.int(nchar(sequence) + 1, size = 3, replace = TRUE))
  deletion_positions <- sort(sample.int(nchar(sequence), size = 3, replace = TRUE))
  
  # Apply mutations, insertions, and deletions
  mutated_sequence <- strsplit(sequence, "")[[1]]
  for (pos in mutation_positions) {
    mutated_sequence[pos] <- sample(nucleotides[nucleotides != mutated_sequence[pos]], size = 1)
  }
  #total_insertions<-0
  for (pos in insertion_positions) {
    insertion_length <- sample(1:4, size = 1)  # Random insertion length (up to 4)
    #total_insertions<-total_insertions+insertion_length
    insertion_chunk <- sample(nucleotides, size = insertion_length, replace = TRUE)
    mutated_sequence <- c(mutated_sequence[1:(pos - 1)], insertion_chunk, mutated_sequence[pos:length(mutated_sequence)])
  }
  #total_deletions<-0
  for (pos in deletion_positions) {
    deletion_length <- sample(1:4, size = 1)  # Random deletion length (up to 4)
    #total_deletions<-total_deletions+deletion_length
    mutated_sequence <- c(mutated_sequence[1:(pos - 1)], mutated_sequence[(pos + deletion_length):length(mutated_sequence)])
  }
  #print(total_insertions)
  #print(total_deletions)
  # Return mutated sequence
  return(paste(mutated_sequence, collapse = ""))
}

V_genes$sequence_new <- lapply(V_genes$sequence,mutate_sequence,mut_avg=mutations)

filepath<-paste0(opt$dir,"/clean.fasta")
file <- file(description = filepath, open = "a")

filepath_naive<-paste0(opt$dir,"/naive.fasta")
file_naive <- file(description = filepath_naive, open = "a")

generate_random_dna_sequence <- function(length) {
  bases <- c("a", "t", "c", "g")
  random_sequence <- paste(sample(bases, length, replace = TRUE), collapse = "")
  return(random_sequence)
}

for (x in 1:10) {
  V<-sample(V_genes$sequence_new,1)
  D<-sample(D_genes$sequence,1)
  J<-sample(J_genes$sequence,1)

  leaves<-opt$leaves
  
  trees<-sim.bd.taxa(n=leaves,numbsim=1,lambda=2.0,mu=0.5)
  tree<-trees[[1]]
  rate<-as.numeric(gsub("_",".",opt$rate))
  junction_rest<-12
  j_half<-floor(junction_rest/2)

 
  junction <- generate_random_dna_sequence(junction_rest)
  
  whole_region <- paste0(tolower(V),substring(junction,1,j_half),tolower(D),substring(junction,j_half+1,junction_rest),tolower(J))
  whole_r_char<-unlist(strsplit(whole_region, ""))

  cat(paste0(">family_",x), file= file_naive, append=TRUE,sep="\n")
  cat(toupper(whole_region), file= file_naive, append=TRUE,sep="\n")

  sims <- simSeq(tree, l=length(whole_r_char), type="DNA",
                  rootseq = whole_r_char, rate=rate)
  
  filepath_fam<-paste0(opt$dir,"/family_",x,".fasta")
  file_fam <- file(description = filepath_fam, open = "a")

  for (s in 1:length(sims)){
    curr_seq<-toupper(paste(as.character(sims[s]), collapse = ''))
    sequence<-curr_seq
    cat(paste0(">family_",x,"_clone_",s), file= file, append=TRUE,sep="\n")
    cat(toupper(sequence), file= file, append=TRUE,sep="\n")
    cat(paste0(">family_",x,"_clone_",s), file= file_fam, append=TRUE,sep="\n")
    cat(toupper(sequence), file= file_fam, append=TRUE,sep="\n")
  }
  close(file_fam)
}
close(file_naive)
close(file)

