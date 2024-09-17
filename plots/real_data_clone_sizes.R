library(data.table)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)

# Define folders and file names
folders <- c("0", "3")
file_name <- "mptp_data_singletons.txt"

# Create an empty list to store data
all_data <- list()

summary_table <- data.frame(TimePoint = character(),
                            Total_Sequences = numeric(),
                            Singleton_Count = numeric(),
                            stringsAsFactors = FALSE)


for (folder in folders) {
  file_path <- paste0("/project/mpennell_978/kavoss/real_data_cf/14007/",folder, "/", file_name)
  
  # Read the TSV file
  data <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Count occurrences of each clone_id
  clone_counts <- data %>%
    group_by(clone_id) %>%
    summarise(count = n())
  
  # Calculate the number of singletons and total sequences
  total_sequences <- nrow(data)
  singleton_count <- sum(clone_counts$count == 1)
  
  # Add to the summary table
  summary_table <- rbind(summary_table, data.frame(TimePoint = folder,
                                                   Total_Sequences = total_sequences,
                                                   Singleton_Count = singleton_count))

  
  # Add the time point to the data
  clone_counts$TimePoint <- folder
  
  # Store the data in the list
  all_data[[folder]] <- clone_counts
}

# Combine all the data into one data frame
combined_data_all <- bind_rows(all_data)

summary_table$percentage_singletons<-summary_table$Singleton_Count/summary_table$Total_Sequences*100
print(summary_table)

m_clones<-ggplot(combined_data_all, aes(x = count, y = TimePoint, color = TimePoint)) +
  geom_jitter(width = 0.2, height = 0.2, size = 1) +
  scale_x_log10() +
  labs(x = "Clone Size (Number of Sequences per Clone)", y = "Time Point", 
       title = "mPTP Clone Size Distribution Across Time Points (log)") +
  theme_minimal() +
  theme(legend.position = "none")


write.csv(combined_data_all, paste0("/project/mpennell_978/kavoss/real_data_cf/14007/","clone_sizes_all.csv"), row.names=FALSE, quote = FALSE)
write.csv(combined_data, paste0("/project/mpennell_978/kavoss/real_data_cf/14007/","clone_sizes.csv"), row.names=FALSE, quote = FALSE)


for (folder in folders) {
  file_path <- paste0("/project/mpennell_978/kavoss/real_data_cf/14007/",folder,"/mapping.csv")
  
  # Read the TSV file
  data <- read.table(file_path, header = TRUE, sep = ",")
  
  # Count occurrences of each clone_id
  clone_counts <- data %>%
    group_by(clone_id) %>%
    summarise(count = n())
  
  # Calculate the number of singletons and total sequences
  total_sequences <- nrow(data)
  singleton_count <- sum(clone_counts$count == 1)
  
  # Add to the summary table
  summary_table <- rbind(summary_table, data.frame(TimePoint = folder,
                                                   Total_Sequences = total_sequences,
                                                   Singleton_Count = singleton_count))
  
  
  # Add the time point to the data
  clone_counts$TimePoint <- folder
  
  # Store the data in the list
  all_data[[folder]] <- clone_counts
}

# Combine all the data into one data frame
combined_data_all <- bind_rows(all_data)

write.csv(combined_data_all, paste0("/project/mpennell_978/kavoss/real_data_cf/14007/","clone_sizes_all_mapping.csv"), row.names=FALSE, quote = FALSE)


vjd<-ggplot(combined_data_all, aes(x = count, y = TimePoint, color = TimePoint)) +
  geom_jitter(width = 0.2, height = 0.2, size = 1) +
  scale_x_log10() +
  labs(x = "Clone Size (Number of Sequences per Clone)", y = "Time Point", 
       title = "V&J+identical CDR3 Clone Size Distribution Across Time Points (log)") +
  theme_minimal() +
  theme(legend.position = "none")

summary_table$percentage_singletons<-summary_table$Singleton_Count/summary_table$Total_Sequences*100

m_clones/vjd

print(summary_table)


