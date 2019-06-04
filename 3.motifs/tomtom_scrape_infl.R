###    General description    ###
### Use this to analyse data following motif detection
### This script gets the useful information out of the tomtom.txt -file. (Basset_motifs.py output) 
### It matches this to TF motif names, the filter influence, and influence per stage
### In conclusion, you'll have one large table with biological
### TF motif names and matching quantitative data
###    End of description    ###

library("plyr")
library("tidyr")

###    Part one: Preparation    ###

# Make sure you have the 4 required files (input and output for Basset_motifs.py) from a given model
# in a single directory, then go there:
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/motif_original_random")
#1. tomtom, this is the main matching table
tomtom <- read.csv("tomtom.txt", header = TRUE, sep = "\t")
colnames(tomtom)[colnames(tomtom) == "X.Query.ID"] <- "Query.ID"

#2. PWMs_adapted.meme, this contains PWM_*number* and biological names
tf_pwms <- readLines("PWMs_adapted.meme")

#3. Filter influence, stdev:
filter_infl <- read.table("table.txt", header = TRUE)

#4. Filter influence per stage:
table_target <- read.table("table_target.txt", header = FALSE)

#Steps to clean and merge merge tomtom and filter influence:
filter_infl$Query.ID <- paste0("filter",row.names(filter_infl))
filter_infl <- filter_infl[, -c(1,3,4)]
tomtom <- tomtom[, c(1,2,6)]
tomtom_infl <- merge(tomtom, filter_infl, by.x = "Query.ID", by.y = "Query.ID", all.x = TRUE)

#Steps to get a biological TF name for each filter
target.ID2motifname <- tf_pwms[grep("MOTIF", tf_pwms)]
target.ID2motifname <- gsub("MOTIF ", "", target.ID2motifname)
target.ID2motifname <- gsub(" \t", "", target.ID2motifname)
target.ID2motifname <- unlist(target.ID2motifname)
target.ID2motifname <- as.data.frame(target.ID2motifname)
target.ID2motifname$target.ID2motifname <- as.character(target.ID2motifname$target.ID2motifname)
target.ID2temp <- strsplit(target.ID2motifname$target.ID2motifname, " ")
target.ID2motifname <- ldply(target.ID2temp)
colnames(target.ID2motifname) <- c("Target.ID", "motifname")
#Merge biological TF names and tomtom_infl.    !!! Warning: you will lose unannotated filters !!!
tomtom_infl_name <- merge(tomtom_infl, target.ID2motifname, by.x = "Target.ID", by.y = "Target.ID", all.x = TRUE)

#This gives you a table of filter influence per stage for all CNN filters
table_target <- table_target[,-2]
table_target <- tidyr::spread(table_target, V3, V4)
colnames(table_target)[1] <- "Query.ID"
table_target$Query.ID <- paste0("filter",table_target$Query.ID)
table_target <- table_target[, c(1,6,3,9,8,7,5,4,2)]

###    End of Part 1    ###

###    Part 2: Final Tables and plots    ###
#Largest table, merge Biological TF names, filter influence and influence per stage
#!!! Warning: you will lose unannoted filter influence per stage rows !!!
ann_complete <- merge(tomtom_infl_name, table_target, by.x = "Query.ID", by.y = "Query.ID", all.x = TRUE)
ann_complete <- ann_complete[, -c(2,4)]
ann_complete$motifname <- gsub(".motif", "" , ann_complete$motifname)

#q-value threshold < 0.05 and std > 0.1 
ann_qselect <- ann_complete[ann_complete$q.value < 0.05, ]
ann_qselect <- ann_qselect[ann_qselect$std > 0.1, ]
rownames(ann_qselect) <- paste(ann_qselect$Query.ID, ann_qselect$motifname, sep = " ")

ann_qselect <- ann_qselect[, -c(1:4)]
ann_qselect <- as.matrix(ann_qselect)
heatmap.2(ann_qselect, col = redblue, dendrogram = "row", Colv = NA, key= T, cexRow = 0.6, margins = c(4,13), density.info = "none", keysize = 1)

#q-value threshold < 0.05 and reduce rows by placing TF names in a single row per filter
ann_qselect2 <- ann_complete[ann_complete$q.value < 0.05, ]
ann_qselect2 <- ann_qselect2[!duplicated(ann_qselect2$Query.ID), ]
ann_qselect2$Query.ID <- gsub("filter", "f", ann_qselect2$Query.ID) 
rownames(ann_qselect2) <- paste(ann_qselect2$Query.ID, ann_qselect2$motifname, sep = " ")
ann_qselect2 <- ann_qselect2[, -c(1:4)]
ann_qselect2 <- as.matrix(ann_qselect2)

heatmap.2(ann_qselect2, col = redblue, dendrogram = "row", Colv = NA, key= T, cexRow = 0.6, margins = c(4,13), density.info = "none", keysize = 1)

###    Make a barplot of std (influence per filter)    ###
bar_select <- ann_complete[, c(1:4)]
bar_select <- bar_select[bar_select$q.value < 0.05, ]
bar_select <- bar_select[!duplicated(bar_select$Query.ID), ]
bar_select$Query.ID <- gsub("filter", "f", bar_select$Query.ID)
bar_select$full_name <- paste(bar_select$Query.ID, bar_select$motifname, sep = " ")
bar_select$full_name <- factor(bar_select$full_name, levels = bar_select$full_name[order(bar_select$std)])
bar_select <- bar_select[, -c(1,2,4)]

p<-ggplot(data=bar_select, aes(x=full_name, y=std, fill = full_name)) +
  geom_bar(stat="identity") + theme_minimal()

# Horizontal bar plot
p + coord_flip()

