#!/usr/bin/env Rscript
library("optparse")
library("Biostrings")
library(data.table)


option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to simulation", metavar="character"),
  make_option(c("-n", "--naive"), type="character", default="out.txt", 
              help="naive.fasta", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

path<-opt$path

## naive seq
fastaFile = readDNAStringSet(opt$naive)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
naive <- data.frame(seq_name, sequence)
print("read in naive file")

for (x in seq_name) {
  discerned<-paste0(path,"/ancestral_sequences/",x,".raxml.ancestralStates.txt")
  if (file.exists(discerned)) {
    print(x)
    disc<-fread(discerned, header=FALSE,col.names=c("node","seq"))
    ancestral_discerned<-tail(disc, n =1)$seq
    file_name<-paste0(path,"/combined/",x,"_combined.fasta")
    file<-file(file_name)
    writeLines(c(">naive",naive[seq_name==x,]$sequence,">PTP_derived",ancestral_discerned), file)
    close(file)
  }
}
print("done")