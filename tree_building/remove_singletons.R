#!/usr/bin/env Rscript
library("optparse")
library("Biostrings")
library(data.table)
require("reshape2")
require("seqRFLP")

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
df <- transform(df, seq_name = colsplit(seq_name, "_", names = c('fam','family', 'c','desc')))
df <- do.call(data.frame, df)
df <- df[,c('seq_name.family','seq_name.desc','sequence')]

## rename columns
colnames(df) <- c("famID", "specific_seqID", "sequence")
counts<-as.data.table(table(df$famID))
singletons<-counts[N==1]$V1

df<-subset(df, !famID %in% singletons)


df$famID <- paste("family", df$famID, "clone",df$specific_seqID, sep="_")
df<-df[,c("famID","sequence")]
dataframe2fas(df, file = paste0(opt$out,"clean", ".fasta"))


