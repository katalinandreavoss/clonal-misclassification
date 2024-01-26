#!/usr/bin/env Rscript
#library("optparse")
library(data.table)
library(tidyverse)
library(patchwork)
library(Biostrings)
library(reshape2)
library(stringr)
library(patchwork)
require(gtools)
clonal_families<-c(4,6,8,10,12,14,16,18, 20)
path<-"/home1/kavoss/simulations/"

get_max <- function(file_path) {
  df<-fread(file_path)
  colnames(df)[1]<-gsub("#taxaorder:","",colnames(df)[1])
  df<-as.data.frame(t(df))
  return(max(df$V1))
}

get_num_fam_no_singletons <- function(file_path) {
  df<-fread(file_path)
  colnames(df)[1]<-gsub("#taxaorder:","",colnames(df)[1])
  df<-as.data.frame(t(df))
  frequencies<-as.data.frame(table(df$V1))
  num<-nrow(frequencies[frequencies$Freq!=1,])
  return(num)
}

get_freq<-function(file_path) {
  df<-fread(file_path)
  colnames(df)[1]<-gsub("#taxaorder:","",colnames(df)[1])
  df<-as.data.frame(t(df))
  frequencies<-as.data.frame(table(df$V1))
  return(frequencies[frequencies$Freq!=1,]$Freq)
}

get_real_family_number<-function(path) {
  clean_data = readDNAStringSet(paste0(path,"clean.fasta"))
  seq_name = names(clean_data)
  sequence = paste(clean_data)
  clean_data<-data.table(seq_name,sequence)
  clean_data<-data.table(colsplit(clean_data$seq_name, "clone", c("family", "clone")),sequence)
  size<-length(unique(clean_data$family))
  return(size)
}



get_median_fam_size<-function(data,df) {
  clone=data$clones
  shm=data$SHM
  leaf=data$leaves
  s=data$sim
  med<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    filter(sim == s)%>%
    select(freq)
  med_real<-median(med$freq)
  return(med_real)
}

get_MSE_median_fam_size<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  real<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    select(median_family_size_real)
  pred<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    select(median_family_size)
  mse<-mean((as.numeric(real$median_family_size_real) - as.numeric(pred$median_family_size))^2)
  return(mse)
}

get_MSE_fam_number<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=data[['leaves']]
  real<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    select(real_fam_number)
  pred<-df%>%
    filter(clones == clone)%>%
    filter(SHM == shm)%>%
    filter(leaves == leaf)%>%
    select(number_families)
  mse<-mean((as.numeric(real$real_fam_number) - as.numeric(pred$number_families))^2)
  return(mse)
}



for (x in clonal_families) {
  print(x)
  curr_path<-paste(path,as.character(x),sep = "")
  filenames <- list.files(path=curr_path, pattern="mega.PTPPartitions.txt", full.names=TRUE,recursive = TRUE)
  
  total_df<-data.frame(filenames)
  total_df$params<-gsub(path,"",total_df$filenames)
  total_df$params<-gsub("mega.PTPPartitions.txt","",total_df$params)
  total_df$params<-str_sub(total_df$params,0,-2)
  
  total_df<-total_df %>%
    separate(params, c("clones","SHM","leaves","sim"), "/")

  total_df$number_families<-as.numeric(lapply(total_df$filenames, get_num_fam_no_singletons))
  
  total_df$freq<-lapply(total_df$filenames, get_freq)
  total_df$filenames<-gsub("mega.PTPPartitions.txt","",total_df$filenames)
  total_df$real_fam_number<-as.numeric(lapply(total_df$filenames, get_real_family_number))
  total_df$median_family_size<-lapply(total_df$freq,median)
  
  
  number_families<-distinct(data.frame(clones=x,SHM=total_df$SHM, leaves=total_df$leaves))
  
  real_freq <-
    filenames <- list.files(path=curr_path, pattern="family_sizes.txt", full.names=TRUE,recursive = TRUE)%>% 
    map_df(~fread(.))
  real_freq$V1<-gsub("_",".",real_freq$V1)
  colnames(real_freq)<-c("clones","SHM","leaves","sim","freq")
  total_df$median_family_size_real<-apply(total_df,1,FUN=get_median_fam_size,df=real_freq)
  number_families$MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=total_df)
  number_families$MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=total_df)
  
  assign(paste("clones",as.character(x),sep = ""),number_families)
  assign(paste("all",as.character(x),sep = ""),total_df)
}

#all_combined<-rbindlist(list(all4,all6, all8,all10,all12,all14,all16,all18,all20))

combined<-rbindlist(list(clones4,clones6, clones8,clones10,clones12,clones14,clones16,clones18,clones20))

#combined$MSE_SHM<-as.numeric(lapply(combined$SHM, get_MSE_number_families, column_name=quo(SHM),df=all_combined))
#combined$MSE_clones<-as.numeric(lapply(combined$clones, get_MSE_number_families, column_name=quo(clones),df=all_combined))
#combined$MSE_leaves<-as.numeric(lapply(combined$leaves, get_MSE_number_families, column_name=quo(leaves),df=all_combined))

#combined$MSE_SHM_family_size<-as.numeric(lapply(combined$SHM, get_MSE_family_size, column_name=quo(SHM),df=all_combined))
#combined$MSE_clones_family_size<-as.numeric(lapply(combined$clones, get_MSE_family_size, column_name=quo(clones),df=all_combined))
#combined$MSE_leaves_family_size<-as.numeric(lapply(combined$leaves, get_MSE_family_size, column_name=quo(leaves),df=all_combined))


combined$clones<-as.character(combined$clones)

MSE_SHM_num_fam<-ggplot(combined, aes(SHM,MSE_num_fam)) + 
  geom_boxplot()

MSE_clones_num_fam<-ggplot(combined, aes(clones,MSE_num_fam)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))


MSE_leaves_num_fam<-ggplot(combined, aes(leaves,MSE_num_fam)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))



MSE_SHM_size_fam<-ggplot(combined, aes(SHM,MSE_fam_size)) + 
  geom_boxplot()

MSE_clones_size_fam<-ggplot(combined, aes(clones,MSE_fam_size)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))

MSE_leaves_size_fam<-ggplot(combined, aes(leaves,MSE_fam_size)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$leaves))


title<-paste(as.character(x), "original families")

MSE_SHM_num_fam+MSE_clones_num_fam+ MSE_leaves_num_fam+plot_annotation('number of discerned families (PTP)',theme=theme(plot.title=element_text(hjust=0.5)))
MSE_SHM_size_fam+MSE_clones_size_fam+ MSE_leaves_size_fam+plot_annotation('median size of discerned families (PTP)',theme=theme(plot.title=element_text(hjust=0.5)))

MSE_SHM_num_fam+MSE_SHM_size_fam+plot_annotation('SHM',theme=theme(plot.title=element_text(hjust=0.5)))
MSE_clones_num_fam+MSE_clones_size_fam+plot_annotation('clones',theme=theme(plot.title=element_text(hjust=0.5)))
MSE_leaves_num_fam+MSE_leaves_size_fam+plot_annotation('clones',theme=theme(plot.title=element_text(hjust=0.5)))



ggplot(combined, aes(clones,MSE_num_fam, fill=leaves)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(combined$clones))


