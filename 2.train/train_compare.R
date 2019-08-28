library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/train_compare/")

results <- list.files(path = " .", pattern = "cat")
for (i in 1:length(results))
  assign(results[i],
         read.csv(results[i], header = FALSE))

compare_input <- read.csv(file = "train_compare.csv", header = TRUE, sep = ',', row.names = 1)
compare_input <- gather(compare_input)
p <- ggplot(compare_input,
            aes(x = key,y = value))
p + geom_boxplot() + theme_minimal() + labs(x = "Input", y = "AUC") +  theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text = element_text(
    angle = 45,
    hjust = 1,
    vjust = 0.9)
)  

# Get the log results in neg_training directory :
#  grep best ./*/log_original.iter* > val_AUCs
# Select best model epoch for each iter

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/neg_train_compare/")
compare_input <- read.csv(file = "val_AUCs", sep = ' ')
