clonal_families<-c(10,20)
SHM<-c("0_001","0_005", "0_01","0_05","0_1","0_2")
leaves <- c("10","20","50","100")
balance <- c("0_0")
sims<-seq(1,50,1)
path<-"/home1/kavoss/panfs/method_comparison"

listed<-tidyr::expand_grid(clonal_families,SHM, leaves,balance,sims)
listed$filenames<-paste(path,listed$clonal_families,listed$SHM,listed$leaves,listed$balance,listed$sims,"analysis_no_singletons.tsv", sep="/")


filenames <- listed$filenames
total_df <- do.call(rbind,lapply(filenames,read.csv,sep="\t"))
  
number_families<-tidyr::expand_grid(unique(total_df$tool),clonal_families,SHM, leaves,balance)
colnames(number_families)<- c("tool","clones","SHM","leaves","balance")
  
number_families$MSE_num_fam<-apply(number_families,1, FUN=get_MSE_fam_number, df=total_df)
number_families$MSE_fam_size<-apply(number_families,1, FUN=get_MSE_median_fam_size, df=total_df)
number_families$MSE_var_fam_size<-apply(number_families,1, FUN=get_MSE_var_fam_size, df=total_df)
number_families$clones<-as.character(number_families$clones)
number_families[number_families$tool=="scoper_hier",]$tool<-"scoper_0_15"
not<-c("mixcr","changeo","scoper_id")
number_families[!number_families$tool %in% not, ]

MSE_SHM_num_fam<-ggplot(number_families[!number_families$tool %in% not, ], aes(SHM,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_SHM_num_fam


MSE_leaves_num_fam<-ggplot(number_families[!number_families$tool %in% not, ], aes(leaves,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(number_families[!number_families$tool %in% not, ]$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_leaves_num_fam

MSE_clones_num_fam<-ggplot(number_families[!number_families$tool %in% not, ]
, aes(clones,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(number_families[!number_families$tool %in% not, ]
$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_clones_num_fam

MSE_num_fam<-ggplot(number_families[!number_families$tool %in% not, ]
                           , aes(tool,MSE_num_fam, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Number of Families (log)")
MSE_num_fam



MSE_SHM_size_fam<-ggplot(number_families[!number_families$tool %in% not, ]
, aes(SHM,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_SHM_size_fam

MSE_leaves_size_fam<-ggplot(number_families[!number_families$tool %in% not, ]
, aes(leaves,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(number_families[!number_families$tool %in% not, ]
$leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_leaves_size_fam

MSE_clones_size_fam<-ggplot(number_families[!number_families$tool %in% not, ]
, aes(clones,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(number_families[!number_families$tool %in% not, ]
$clones))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_clones_size_fam

MSE_size_fam<-ggplot(number_families[!number_families$tool %in% not, ]
                            , aes(tool,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")
MSE_size_fam

MSE_clones_num_fam+ MSE_leaves_num_fam+plot_annotation('number of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")
MSE_clones_size_fam+ MSE_leaves_size_fam+plot_annotation('median size of discerned families',theme=theme(plot.title=element_text(hjust=0.5)))+ plot_layout(guides = "collect")




test<-number_families[!number_families$tool %in% not, ]


ggplot(test[test$SHM=="0_001", ]
                            , aes(leaves,MSE_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_x_discrete(limits = unique(number_families[!number_families$tool %in% not, ]
                                   $leaves))+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Median Family Size (log)")


true<-total_df[c("clones" ,  "SHM", "leaves", "balance","real_med_size")]
colnames(true)<-c("clones" ,  "SHM", "leaves", "balance","median_size_families")
true$tool<-"true"

combined_data<-total_df[c("tool","clones" ,  "SHM", "leaves", "balance","median_size_families")]
only<-c("true","scoper_hier")
combined_data<-rbind(combined_data,true)  
combined_data$leaves<-as.character(combined_data$leaves)

ggplot(combined_data[combined_data$tool %in% only, ]
       , aes(leaves,median_size_families, fill=tool)) + 
  geom_boxplot()+
  facet_grid(SHM ~ .)+
  scale_y_continuous(trans='log10')+
  ylab("Median Family Size (log)")
  
  
  
  
  
