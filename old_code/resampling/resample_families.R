#!/usr/bin/env Rscript
## resample the sequences in a family at a certain rate
## this is for one simulation or for one real life data set

## load packages
require("Biostrings")
require("reshape2")
library("devtools")
#install_github("helixcn/seqRFLP")
require("seqRFLP")
library("optparse")

option_list = list(
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="path to fasta file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="path to output directory", metavar="character"),
  make_option(c("-s", "--shuffle"), type="character", 
              help="rate of shuffle",metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path<-opt$fasta
#opt$out <-"/home1/kavoss/simulation_from_scratch/sim_10/0_01/shuffle/"
#path <-"/home1/kavoss/simulation_from_scratch/sim_10/0_01/all.fasta"

## load fasta file
fastaFile <- readDNAStringSet(path)
seq_name <- names(fastaFile)
sequence <- paste(fastaFile)

## convert to data frame
df <- data.frame(seq_name, sequence)

## separate columns for seqID, famID, seqID within family
df <- transform(df, seq_name = colsplit(seq_name, "_", names = c('fam','family', 'c','desc')))

df <- do.call(data.frame, df)
df <- df[,c('seq_name.family','seq_name.desc','sequence')]

## rename columns
colnames(df) <- c("famID", "specific_seqID", "sequence")

## set r to the desired resampling level
shuffle <- gsub("\\_", "\\.", opt$shuffle)
r <- as.numeric(shuffle)

## sample value 1 between 0 and 1
u <- list()
for (i in 1:dim(df)[[1]]){
  u[[i]] <- runif(1)
}
df$u <- u

## if u <= r then move the sequence to a new family
df$fam_new <- 0

for (i in 1:dim(df)[[1]]){
  if (df$u[i] <= r) {
    fam_new <- as.vector(sample(df$famID))[1]
    df$fam_new[i] <- fam_new
    }
  else{
    fam_old <- as.vector(df$famID)[i]
    df$fam_new[i] <- fam_old
  }
}

## make new dataframes for each new family
dfs <- split(df, f = df$fam_new)
dfs2 <- list()
dfs3 <- list()

## extract only information about old and new family + sequence
for (i in 1:length(dfs)){
df2 <- data.frame(lapply(dfs[[i]], as.character), stringsAsFactors=FALSE)
dfs2[[i]] <- df2
dfs2[[i]]$id <- paste("family", dfs2[[i]]$famID, "clone",dfs2[[i]]$specific_seqID, "new_family",dfs2[[i]]$fam_new, sep="_")
df3 <- data.frame(dfs2[[i]]$id, dfs2[[i]]$sequence)
dfs3[[i]] <- df3
}

## save files, one fasta file for each family
lapply(1:length(dfs3), function(i) dataframe2fas(dfs3[[i]], 
                                                 file = paste0(opt$out,'/family_',names(dfs)[i], ".fasta")))



