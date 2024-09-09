#!/usr/bin/env Rscript
library("optparse")
library("Biostrings")
library(data.table)
require("reshape2")
require("phangorn")
require("seqRFLP")
library(stringdist)
library(parallel)


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

# Number of cores to use
num_cores <- parallel::detectCores() - 1


# Split the sequences into chunks to distribute across cores
sequence_chunks <- split(sequences, rep(1:num_cores, length.out = length(sequences)))

pairwise_levenshtein <- function(seq_subset, all_sequences, chunk_num) {
  message("Processing chunk ", chunk_num, " with ", length(seq_subset), " sequences...")
  
  # Initialize an empty matrix to store distances
  dist_matrix_chunk <- matrix(0, nrow = length(seq_subset), ncol = length(all_sequences))
  
  # Fill in the distance matrix
  for (i in seq_along(seq_subset)) {
    dist_matrix_chunk[i, ] <- stringdist::stringdist(seq_subset[i], all_sequences, method = "lv")
  }
  
  message("Finished chunk ", chunk_num)
  return(dist_matrix_chunk)
}

dist_list <- mclapply(1:length(sequence_chunks), function(i) {
  pairwise_levenshtein(sequence_chunks[[i]], sequences, i)
}, mc.cores = num_cores)

# Combine the results into a single matrix
dist_matrix <- do.call(rbind, dist_list)


#dist_matrix <- stringdist::stringdistmatrix(sequences, method = "lv")
dist_matrix <- as.matrix(dist_matrix)
rownames(dist_matrix) <- df$seq_name
colnames(dist_matrix) <- df$seq_name

tree <- upgma(dist_matrix)

write.tree(tree, file = paste0(opt$out,"/upgma.tre"))





