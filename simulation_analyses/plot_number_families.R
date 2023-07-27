library(data.table)
library(tidyverse)
library(patchwork)


clonal_families<-c(2,4,6,8)

path<-"/home1/kavoss/simulation_from_scratch/sim_"

get_max <- function(file_path) {
  df<-fread(file_path)
  colnames(df)[1]<-gsub("#taxaorder:","",colnames(df)[1])
  df<-as.data.frame(t(df))
  return(max(df$V1))
}

get_freq<-function(file_path) {
  df<-fread(file_path)
  colnames(df)[1]<-gsub("#taxaorder:","",colnames(df)[1])
  df<-as.data.frame(t(df))
  frequencies<-as.data.frame(table(df$V1))
  return(frequencies$Freq)
}

for (x in clonal_families) {
  print(x)
  curr_path<-paste(path,as.character(x),sep = "")
  filenames <- list.files(path=curr_path, pattern="all_PTP.PTPPartitions.txt", full.names=TRUE,recursive = TRUE)
  total_df<-data.frame(filenames)
  SHM<-gsub(curr_path,"",total_df$filenames)
  SHM<-gsub("all_PTP.PTPPartitions.txt","",SHM)
  SHM<-gsub("/","",SHM)
  SHM<-gsub("_",".",SHM)
  total_df$SHM<-SHM
  
  total_df$number_families<-as.numeric(lapply(total_df$filenames, get_max))
  total_df$freq<-lapply(total_df$filenames, get_freq)
  
  freq_table<-total_df[,c("SHM","freq")]
  freq_table<-freq_table %>% unnest(freq)
  
  title<-paste(as.character(x), "original families")
  p<-ggplot(total_df, aes(SHM,number_families)) + 
    geom_point()+
    geom_hline(yintercept=x, linetype="dashed", color = "red")+
    ylim(0, 450)+
    ggtitle(title)

  m<-ggplot(freq_table, aes(SHM,freq)) + 
    geom_boxplot()+
    ylim(0, 200)+
    ggtitle(title)
  
  assign(paste("p",as.character(x),sep = ""),p)
  assign(paste("m",as.character(x),sep = ""),m)
}


(p2+p4)/(p6+p8)+ plot_annotation('number of discerned families (PTP)',theme=theme(plot.title=element_text(hjust=0.5)))

(m2+m4)/(m6+m8)+ plot_annotation('size of discerned families (PTP)',theme=theme(plot.title=element_text(hjust=0.5)))




