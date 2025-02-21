library(dplyr)
library(tidyverse)
library(clue)  # For Adjusted Rand Index

# Load both files
path<-"/project/mpennell_978/kavoss/real_data_cf/14007/"
folder<-"0"
mapping_file = paste0("mapping.csv")
singleton_file = "mptp_data_singletons.txt"

mapping_file <- read.csv(paste0(path,folder,"/",mapping_file))
singleton_file <- read.delim(paste0(path,folder,"/",singleton_file), sep = "\t")


mapping_group_sizes <- mapping_file %>%
  group_by(clone_id) %>%
  summarise(size = n())%>%
  mutate(dataset = "VJCDR3")

singleton_group_sizes <- singleton_file %>%
  group_by(clone_id) %>%
  summarise(size = n())%>%
  mutate(dataset = "mPTP")

merged <- bind_rows(mapping_group_sizes, singleton_group_sizes)

ggplot(merged, aes(x = dataset, y = size, fill = dataset)) +
  geom_boxplot() +
  labs(title = "Comparison of Clone Sizes Between Datasets",
       x = "Dataset",
       y = "Clone Size (Number of Sequences)") +
  scale_y_log10() +  # Log scale for better visibility of differences
  theme_minimal()


f1<-read.csv(paste0(path,folder,"/sensitivity_precision.tsv"), sep = "\t")

ggplot(f1, aes( y=rand))+
  geom_boxplot()+
  ylab('Rand Index')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )

ggplot(f1, aes( y=f1))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )


ggplot(f1, aes( y=sensitivity))+
  geom_boxplot()+
  ylab('F1-Score')+
  scale_fill_brewer(palette = "Paired")+
  theme(text = element_text(face = "bold",size = 20),
        panel.background = element_rect(fill = "white") )


### scatter plot

mapping_df <- mapping_file 
colnames(mapping_df)<-c("sequence_id","family")


df <- merge(singleton_file, mapping_df, by = "sequence_id", all.x = TRUE)

clone_counts <- df %>%
  group_by(clone_id) %>%
  summarise(clone_count = n()) 

family_counts <- df %>%
  group_by(family) %>%
  summarise(family_count = n()) 

colnames(family_counts)[2]<-"family_id_count"


# Merge counts with the original data
df <- df %>%
  left_join(clone_counts, by = "clone_id") %>%
  left_join(family_counts, by = "family")

write.csv(df, paste0(path,folder,"/grouping_merged"), row.names=FALSE, quote = FALSE)


# Plot the data
ggplot(df, aes(x = family_id_count, y = clone_count)) +
  geom_point(alpha = 0.1,size=0.5) +
  labs(x = "Number of Sequences with Same V&J group and identical CDR3",
       y = "Number of Sequences with Same mPTP group") +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal()

ggplot(df, aes(x = family_id_count, y = clone_count)) +
  geom_point(alpha = 0.5,size=0.5) +
  labs(x = "Number of Sequences with Same V&J group and identical CDR3",
       y = "Number of Sequences with Same mPTP group") +
  theme_minimal()

