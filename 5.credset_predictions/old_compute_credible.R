###   Quick Description   ###

#THIS SCRIPT IS TOO MESSY, BUT CONTAINS USEFULL SNIPPETS.
#This is an old script, used to analyse fasta files of variants, their scores and PPa values.
#This first script takes "HRC_all_variants.fasta" and "HRC_credset.snp_ann.txt" as input
#The former contains names of variants and their fasta seqs in one column
#The latter contains names and locations in the hg19 genome


HRC_all_variants <- read.table("HRC.all_variants.fasta")
name_to_loc <- read.table("HRC_credset.snp_ann.txt")

#Merge these and then merge with results from "credible_set_better_sep.txt", which contains the CNN-calculated likelihood
#for each of the fasta seqs rows (in the same order!)
nuc_seqs <- HRC_all_variants[seq(2,nrow(HRC_all_variants),2), ]
hrc_all_varnames <- HRC_all_variants[seq(1,nrow(HRC_all_variants),2), ]
hrc_all_name_seq <- data.frame(hrc_all_varnames, nuc_seqs)
colnames(hrc_all_name_seq) <- c("name_and_var", "sequences")

#Duplicate name column for merge
hrc_all_name_seq$name_only <- hrc_all_name_seq$name_and_var
hrc_all_name_seq <- hrc_all_name_seq[ , c(3,1,2)]
hrc_all_name_seq$name_only <- gsub(">", "", hrc_all_name_seq$name_only)

#Probably not the most efficient way:
hrc_all_name_seq$name_only <- gsub("_ref", "", hrc_all_name_seq$name_only)
hrc_all_name_seq$name_only <- gsub("_alt", "", hrc_all_name_seq$name_only)

#Tailor name_to_loc
name_to_loc2 <- name_to_loc[, -3]
colnames(name_to_loc2) <- c("name_only", "gen_loc")

# CNN from credible_set_better_sep.txt
cred_set_results <- read.table("credible_set_sep.txt")

### Make new dataframe for CNN prediction, this time only unique names and sequences ###
#make it like hrc_all_variants, but unique.
name_and_seq <- hrc_all_name_seq[, -1]
uniques <- name_and_seq[!duplicated(name_and_seq), ]
library(dplyr)
to_cnn <- as.vector(t(uniques))
to_cnn <- as.data.frame(to_cnn)
write.csv(to_cnn, file="HRC.all_variants_new.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

### Redo this part when you finish training ###
#You should determine the score differences and merge that dataframe with the one containing locations

### Try this while waiting for GPU ###
#load libraries
library(gplots)
library(RColorBrewer)
#Select first 6000 names from name and sequence list and from cred_set_results
temp_selected_names <-hrc_all_name_seq[1:6000, 1:2 ]
temp_selected_CNN <- cred_set_results[1:6000]
#Calculate differences from cred_set_results (subtraction) and make plot
temp_diff <- temp_selected_CNN[seq(1,nrow(temp_selected_CNN),2), ] - temp_selected_CNN[seq(2,nrow(temp_selected_CNN),2), ]
temp_diff <- as.matrix(temp_diff)
pdf('heat.pdf')
heatmap.2(temp_diff, col = brewer.pal(9,"Blues"), margins = c(14, 14), density.info = "none", dendrogram = "none", Colv ="FALSE", trace= "none")
dev.off()
#Concatenate name and difference
single_names <- as.data.frame(temp_selected_names$name_only[seq(1,nrow(temp_selected_names),2)])
name_and_dif <- cbind(single_names,temp_diff)
#Find gene location in name_to_loc



### It goes downhill here ###
#DOESN'T MERGE IN ORIGINAL ORDER
name_loc_seq <- merge(name_to_loc2, hrc_all_name_seq)
#TRY TO MERGE IN ORDER SO YOU CAN COMPARE TO CNN DIFFERENCES

#Other problem: the results from CNN and the fasta input don't have the same dimensions:
#> dim(cred_set_results)
#[1] 219558      8
#> dim(nuc_seqs)
#[1] 219598
### That's a 40 seq difference! ###
### Remove duplicates and redo Basset CNN prediction ###
check <- hrc_all_name_seq[duplicated(hrc_all_name_seq), ]
name_and_seq <- hrc_all_name_seq[, -1]
uniques <- name_and_seq[!duplicated(name_and_seq), ]
uniques_new <- uniques[rep(1:nrow(uniques), 1, each=2), ]
seq(2,nrow(uniques_new),2)

#Some additional remarks:
#There are duplicates in names of sequences (hrc_all_varnames) and in the sequences themselves
#There are 40 duplicate names, and 538 duplicate sequences. None of the duplicates are var and ref of the same name.
