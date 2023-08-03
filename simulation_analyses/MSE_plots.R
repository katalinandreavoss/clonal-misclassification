#!/usr/bin/env Rscript
#library("optparse")
library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)

clonal_families<-c(6,8,10,12,14,16,18, 20)
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

get_real_family_size<-function(path) {
  clean_data = readDNAStringSet(paste0(path,"clean.fasta"))
  seq_name = names(clean_data)
  sequence = paste(clean_data)
  clean_data<-data.table(seq_name,sequence)
  clean_data<-data.table(colsplit(clean_data$seq_name, "clone", c("family", "clone")),sequence)
  size<-length(unique(clean_data$family))
  return(size)
}

get_MSE_number_families<-function(value, column_name, df) {
  pred<-as.numeric(unlist(df%>%
    filter(!!column_name ==value)%>%
    select(number_families)))
  real<-as.numeric(unlist(df%>%
    filter(!!column_name ==value)%>%
    select(real_fam_size)))
  mse<-mean((real - pred)^2)
  return(mse)
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
  shm<-param_table$SHM
  leaves<-param_table$leaves
  
  total_df$SHM<-shm
  total_df$leaves<-leaves
  total_df$number_families<-as.numeric(lapply(total_df$filenames, get_max))
  total_df$freq<-lapply(total_df$filenames, get_freq)
  total_df$filenames<-gsub("mega.PTPPartitions.txt","",total_df$filenames)
  total_df$real_fam_size<-as.numeric(lapply(total_df$filenames, get_real_family_size))
  total_df$clones<-x
  
  number_families<-distinct(data.frame(clones=x,SHM=total_df$SHM, leaves=total_df$leaves))
  
  assign(paste("clones",as.character(x),sep = ""),number_families)
  assign(paste("all",as.character(x),sep = ""),total_df)
}

all_combined<-rbindlist(list(all6, all8,all10,all12,all14,all16,all18,all20))

combined<-rbindlist(list(clones6, clones8,clones10,clones12,clones14,clones16,clones18,all20))

combined$MSE_SHM<-as.numeric(lapply(combined$SHM, get_MSE_number_families, column_name=quo(SHM),df=all_combined))
combined$MSE_clones<-as.numeric(lapply(combined$clones, get_MSE_number_families, column_name=quo(clones),df=all_combined))
combined$MSE_leaves<-as.numeric(lapply(combined$leaves, get_MSE_number_families, column_name=quo(leaves),df=all_combined))

combined$clones<-as.character(combined$clones)

ggplot(combined, aes(SHM,MSE_SHM)) + 
  geom_boxplot()

ggplot(combined, aes(clones,MSE_clones)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))

ggplot(combined, aes(leaves,MSE_leaves)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))

freq_table<-total_df[,c("SHM","freq","leaves")]
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



m<-ggplot(freq_total, aes(SHM,freq, fill=data)) + 
  geom_boxplot()+
  ggtitle(title)

assign(paste("p",as.character(x),sep = ""),p)
assign(paste("m",as.character(x),sep = ""),m)



