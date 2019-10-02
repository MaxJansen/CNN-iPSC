# Convert the extracted test set from hdf5 file to a proper test set act
# It was an array with a strange numerical format. This is necessary because 
# HDF5 permute selects random peaks for test set creation
library(stringr)

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/better_train/prepare_test/")

test_set_act <- read.table("test_act.txt")
example_act <- read.table("original_islets_act.txt")
sample_names <- read.table("8col1CPM_samples.txt")
bed_coord <- read.table("test_coord.txt")


# Make the act-file that will be used to assess prediction
bed_coord$V1 <- sub(pattern = 'b', replacement = '', x = bed_coord$V1)
bed_coord$V1 <- gsub("'", '', bed_coord$V1)

row.names(test_set_act) <- bed_coord$V1
colnames(test_set_act) <- sample_names$V1

write.table(x = test_set_act, file = "final_test_set_act.txt", quote = FALSE, sep = "\t", col.names = TRUE)

# Make the bed file that will be used to get fasta
new_cols <- as.data.frame(str_split_fixed(bed_coord$V1, ":", 2))
new_coord <- as.data.frame(str_split_fixed(new_cols$V2, "-", 2))
new_coord$chrom <- new_cols$V1
new_coord$name <- "."
new_coord$score <- 1
new_coord$strand <- "+"
new_coord$V2 <- gsub("\\(\\+)", '', new_coord$V2)
new_coord <- new_coord[, c(3,1,2,4,5,6)]

write.table(x = new_coord, file = "final_test_set_bed.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)         