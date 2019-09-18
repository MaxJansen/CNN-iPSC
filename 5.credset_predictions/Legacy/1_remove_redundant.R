###   Quick Description   ###

#This first script takes "HRC_all_variants.fasta" as input
#The former contains names of variants and their fasta seqs in one column
#It removes redundancy and creates a new table, which can be used for predictions by a neural network

setwd("~/Oxford 2.0/Scripts/CNN_project/Data/credset_predictions/")
HRC_all_variants <- read.table("HRC.all_variants_new.fasta")

#Split HRC_all_variants into names and sequences
nuc_seqs <- HRC_all_variants[seq(2,nrow(HRC_all_variants),2), ]
hrc_all_varnames <- HRC_all_variants[seq(1,nrow(HRC_all_variants),2), ]
hrc_all_name_seq <- data.frame(hrc_all_varnames, nuc_seqs)
colnames(hrc_all_name_seq) <- c("name_and_var", "sequences")

#Duplicate name column for merge in later steps
hrc_all_name_seq$name_only <- hrc_all_name_seq$name_and_var
hrc_all_name_seq <- hrc_all_name_seq[ , c(3,1,2)]
hrc_all_name_seq$name_only <- gsub(">", "", hrc_all_name_seq$name_only)
hrc_all_name_seq$name_only <- gsub("_ref", "", hrc_all_name_seq$name_only)
hrc_all_name_seq$name_only <- gsub("_alt", "", hrc_all_name_seq$name_only)

### Make new dataframe for CNN prediction, this time only unique names and sequences ###
#make it like hrc_all_variants, but unique.
name_and_seq <- hrc_all_name_seq[, -1]
uniques <- name_and_seq[!duplicated(name_and_seq), ]
to_cnn <- as.vector(t(uniques))
to_cnn <- as.data.frame(to_cnn)
write.csv(to_cnn, file="HRC.all_variants_new.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Make new dataframe with unique names, names and vars and sequences. This is useful for merging later!
unique_all_name_seq <- hrc_all_name_seq[!duplicated(hrc_all_name_seq), ]
write.csv(unique_all_name_seq, file="unique_all_name_seq.csv", row.names = FALSE, quote = FALSE)

#Quick test
test <- read.table("unique_all_name_seq.csv", header = TRUE, sep = ",")
