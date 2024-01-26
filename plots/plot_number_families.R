#!/usr/bin/env Rscript
#library("optparse")
library(data.table)
library(tidyverse)
library(patchwork)

#option_list = list(
#  make_option(c("-p", "--path"), type="character", default=NULL, 
#              help="path to simulations", metavar="character"),
#  make_option(c("-o", "--out"), type="character", default="out.txt", 
#              help="output file name [default= %default]", metavar="character")
#); 

#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);


clonal_families<-c(8,6,14)
#opt$path
#path<-paste(opt$path,"sim_",sep = "")
path<-"/home1/kavoss/simulations/"


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
  filenames <- list.files(path=curr_path, pattern="mega.PTPPartitions.txt", full.names=TRUE,recursive = TRUE)
  total_df<-data.frame(filenames)
  params<-gsub(curr_path,"",total_df$filenames)
  params<-gsub("mega.PTPPartitions.txt","",params)
  params<-data.table(strsplit(params,split="/",fixed=T))
  
  param_table <- params %>% 
    cbind(., do.call('rbind', .$V1)) %>% 
    select(-V1)
  colnames(param_table)<-c("SHM","leaves","sim")
  SHM<-param_table$SHM
  SHM<-gsub("_",".",SHM)
  leaves<-param_table$leaves
  
  total_df$SHM<-SHM
  total_df$leaves<-leaves
  total_df$number_families<-as.numeric(lapply(total_df$filenames, get_max))
  total_df$freq<-lapply(total_df$filenames, get_freq)
  
  freq_table<-total_df[,c("SHM","freq")]
  freq_table<-freq_table %>% unnest(freq)

  real_freq <-
    filenames <- list.files(path=curr_path, pattern="family_sizes.txt", full.names=TRUE,recursive = TRUE)%>% 
    map_df(~fread(.))
  real_freq$V1<-gsub("_",".",real_freq$V1)
  colnames(real_freq)<-c("clones","SHM","leaves","sim","freq")
  freq_total<-rbind(freq_table,real_freq[,c("SHM","freq")])
  freq_total<-rbindlist(list(freq_table, real_freq[,c("SHM","freq")]), idcol = 'data')[, data:= paste0('set', data)][]
  freq_total[data=="set1"]$data<-"discerned families"
  freq_total[data=="set2"]$data<-"real families"
  
  title<-paste(as.character(x), "original families")
  p<-ggplot(total_df, aes(SHM,number_families)) + 
    geom_point()+
    geom_hline(yintercept=x, linetype="dashed", color = "red")+
    ggtitle(title)

  m<-ggplot(freq_total, aes(SHM,freq, fill=data)) + 
    geom_boxplot()+
    ggtitle(title)
  
  assign(paste("p",as.character(x),sep = ""),p)
  assign(paste("m",as.character(x),sep = ""),m)
}


number<-(p6)/(p8+p14)+ plot_annotation('number of discerned families (PTP)',theme=theme(plot.title=element_text(hjust=0.5)))

size<-(m6)/(m8+m14)+ plot_annotation('size of discerned families (PTP)',theme=theme(plot.title=element_text(hjust=0.5)))
path<-"/home1/kavoss/simulation_from_scratch/"
ggsave(paste(path,'numer_disc_families.png',sep = ""), number)
ggsave(paste(path,'size_disc_families.png',sep = ""), size)


