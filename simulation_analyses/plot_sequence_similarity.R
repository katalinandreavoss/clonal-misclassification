library(data.table)
library(tidyverse)
library(patchwork)
library("Biostrings")
library(RecordLinkage)
library(ggplot2)
library(patchwork)

clonal_families<-c(2,4,6,8)

path<-"/home1/kavoss/simulation_from_scratch/sim_"


get_similarity<-function(file_path) {
  fastaFile = readDNAStringSet(file_path)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  similarity<-levenshteinSim(sequence[1], sequence[2])
  return(similarity*100)
}

for (x in clonal_families) {
  print(x)
  curr_path<-paste(path,as.character(x),sep = "")
  filenames <- list.files(path=curr_path, pattern="combined_aligned.fasta", full.names=TRUE,recursive = TRUE)
  total_df<-data.frame(filenames)
  total_df$similarity<-as.numeric(lapply(total_df$filenames, get_similarity))
  title<-paste(as.character(x), "original families")
  p<-ggplot(total_df, aes(similarity)) + 
    geom_histogram(binwidth=1)+
    ggtitle(title)+
    coord_cartesian(xlim = c(0, 100))
  assign(paste("p",as.character(x),sep = ""),p)
}

(p2+p4)/(p6+p8)+ plot_annotation('percent sequence similarity compared to true ancestral sequence',theme=theme(plot.title=element_text(hjust=0.5)))


