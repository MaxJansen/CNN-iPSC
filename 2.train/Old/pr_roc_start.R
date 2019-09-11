library(precrec)
library(ggplot2)
library(grid)
library(tidyr)
library(plyr)
library(dplyr)

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_train/iter1")

temp = list.files(pattern="roc*")
for (i in 1:length(temp)) assign(temp[i], read.delim(temp[i], header = FALSE))
selected <- ls()[grep('roc', ls())]


roc1.txt$name <- "BLC"
roc2.txt$name <- "DE"
roc3.txt$name <- "EN"
roc4.txt$name <- "EP"
roc5.txt$name <- "GT"
roc6.txt$name <- "PE"
roc7.txt$name <- "PF"
roc8.txt$name <- "iPSC"

roc_all.txt <- dplyr::bind_rows(roc1.txt, roc2.txt, roc3.txt, roc4.txt, roc5.txt, roc6.txt, roc7.txt, roc8.txt)

ggplot(data=roc_all.txt, aes(x=V1, y=V2, group=name, color=name)) +
  geom_line() + 
  scale_color_brewer(palette="Paired")+
  labs(title="ROC plots of model test set",x="FPR", y = "TPR") +
  theme_minimal()