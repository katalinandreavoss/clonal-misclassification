combined<-read.csv("/home1/kavoss/panfs/method_comparison/output_without_singletons.csv")
all_combined<-read.csv("/home1/kavoss/panfs/method_comparison/output_all_without_singletons.csv")

get_MSE_var_fam_size<-function(data,df) {
  clone=data[['clones']]
  shm=data[['SHM']]
  leaf=as.numeric(data[['leaves']])
  balances=data[['balance']]
  which_tool=data[['tool']]
  real<-df%>%
    filter(clones == str_trim(clone))%>%
    filter(SHM == str_trim(shm))%>%
    filter(leaves == str_trim(leaf))%>%
    filter(balance == str_trim(balances))%>%
    filter(tool == str_trim(which_tool))%>%
    select(real_var_size)
  pred<-df%>%
    filter(clones == str_trim(clone))%>%
    filter(SHM == str_trim(shm))%>%
    filter(leaves == str_trim(leaf))%>%
    filter(balance == str_trim(balances))%>%
    filter(tool == str_trim(which_tool))%>%
    select(var_size_families)
  vector<-as.numeric(real[!is.na(pred)]) - as.numeric(pred[!is.na(pred)])
  mse<-mean(vector[!is.nan(vector)]^2)
  return(mse)
}

combined$MSE_var_fam_size<-apply(combined,1, FUN=get_MSE_var_fam_size, df=all_combined)

combined$clones<-as.character(combined$clones)
combined$leaves<-as.character(combined$leaves)
combined$leaves <- factor(combined$leaves, levels=c("10","20","50","100"))

combined$population_size<-as.numeric(combined$clones)*as.numeric(as.character(combined$leaves))

MSE_SHM_var_fam<-ggplot(combined, aes(SHM,MSE_var_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Variance Family Size (log)")
MSE_SHM_var_fam

MSE_clones_var_fam<-ggplot(combined, aes(clones,MSE_var_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Variance Family Size (log)")
MSE_clones_var_fam

MSE_leaves_var_fam<-ggplot(combined, aes(leaves,MSE_var_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Variance Family Size (log)")
MSE_leaves_var_fam

MSE_bal_var_fam<-ggplot(combined, aes(balance,MSE_var_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Variance Family Size (log)")
MSE_bal_var_fam

MSE_var_fam<-ggplot(combined, aes(tool,MSE_var_fam_size, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("Mean Squared Error - Variance Family Size (log)")
MSE_var_fam



real_data<-all_combined[c("clones","SHM","leaves","balance","real_num_fam","real_med_size","real_var_size")]
real_data<-unique(real_data)
real_data$tool<-"real"
colnames(real_data)<-c("clones","SHM","leaves","balance","num_families","median_size_families","var_size_families","tool")

pred_data<-all_combined[c("tool","clones","SHM","leaves","balance","num_families","median_size_families","var_size_families")]

all_data<-rbind(real_data,pred_data)
write.csv(all_data, "/home1/kavoss/panfs/method_comparison/output_melted.csv", row.names=FALSE)

qqPlot(all_data,)



all_combined$clones<-as.character(all_combined$clones)
all_combined$leaves<-as.character(all_combined$leaves)
all_combined$leaves <- factor(all_combined$leaves, levels=c("10","20","50","100"))


singletons_clones<-ggplot(all_combined, aes(clones,singletons, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("#singletons (log)")
singletons_clones

singletons_SHM<-ggplot(all_combined, aes(SHM,singletons, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("#singletons (log)")
singletons_SHM

singletons_leaves<-ggplot(all_combined, aes(leaves,singletons, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("#singletons (log)")
singletons_leaves

singletons_balance<-ggplot(all_combined, aes(balance,singletons, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("#singletons (log)")
singletons_balance

singletons<-ggplot(all_combined, aes(tool,singletons, fill=tool)) + 
  geom_boxplot()+
  scale_y_continuous(trans='log10')+
  ylab("#singletons (log)")
singletons



