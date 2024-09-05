#!/usr/bin/env Rscript
library("optparse")
library("Biostrings")
library(data.table)
require("reshape2")
require("phangorn")
require("seqRFLP")
library(stringdist)


option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="path to simulation", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

fastaFile = readDNAStringSet(opt$fasta)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)

sequences <- df$sequence

# Calculate pairwise Levenshtein distances
dist_matrix <- stringdist::stringdistmatrix(sequences, method = "lv")
dist_matrix <- as.matrix(dist_matrix)
rownames(dist_matrix) <- df$seq_name
colnames(dist_matrix) <- df$seq_name

tree <- upgma(dist_matrix)

write.tree(tree, file = paste0(opt$out,"/upgma.tre"))





