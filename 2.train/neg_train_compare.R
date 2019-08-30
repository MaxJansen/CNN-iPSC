library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/neg_train_compare/")

compare_input <- read.csv(file = "neg_train_compare.csv", header = TRUE, sep = ',', row.names = 1)
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


