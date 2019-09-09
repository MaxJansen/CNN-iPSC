###   Quick Description   ###
#After compute_credible to preprocess the tables, this script can be used to
#make graphs and basic statistics

library(dplyr)
setwd("~/Oxford/RealScripts/credible_sets/data")

cred_set_results <- read.table("credible_set_sep.txt")
name_var_seq <- read.table("unique_all_name_seq.csv", header = TRUE, sep = ",")
name_to_loc <- read.table("HRC_credset.snp_ann.txt")


#Calculate the differences for each locus and each stage
#You take alternating rows and subtract the second from the first.
#The first row is the 
diff <- cred_set_results[seq(1,nrow(cred_set_results),2), ] - cred_set_results[seq(2,nrow(cred_set_results),2), ]


#Take the "diff" values of the columns(stages) and make a scatterplot for each vs. each
plot(diff$V1, diff$V2, asp = 1)
plot(diff$V1, diff$V3, asp = 1)
plot(diff$V1, diff$V4, asp = 1)
plot(diff$V1, diff$V5, asp = 1)
plot(diff$V1, diff$V6, asp = 1)
plot(diff$V1, diff$V7, asp = 1)
plot(diff$V1, diff$V8, asp = 1)

plot(diff$V2, diff$V3, asp = 1)
plot(diff$V2, diff$V4, asp = 1)
plot(diff$V2, diff$V5, asp = 1)
plot(diff$V2, diff$V6, asp = 1)
plot(diff$V2, diff$V7, asp = 1)
plot(diff$V2, diff$V8, asp = 1)

plot(diff$V3, diff$V4, asp = 1)
plot(diff$V3, diff$V5, asp = 1)
plot(diff$V3, diff$V6, asp = 1)
plot(diff$V3, diff$V7, asp = 1)
plot(diff$V3, diff$V8, asp = 1)

plot(diff$V4, diff$V5, asp = 1)
plot(diff$V4, diff$V6, asp = 1)
plot(diff$V4, diff$V7, asp = 1)
plot(diff$V4, diff$V8, asp = 1)

plot(diff$V5, diff$V6, asp = 1)
plot(diff$V5, diff$V7, asp = 1)
plot(diff$V5, diff$V8, asp = 1)

plot(diff$V6, diff$V7, asp = 1)
plot(diff$V6, diff$V8, asp = 1)

plot(diff$V7, diff$V8, asp = 1)


#To compare diffs to position PPA, merge dataframes
#(Optional) Tailor name_to_loc
name_to_loc2 <- name_to_loc[, -3]
colnames(name_to_loc2) <- c("name_only", "gen_loc")

#Add names and locations to diff
diff_name <- diff
diff_name$name <- name_var_seq$name_only[seq(1,nrow(name_var_seq),2)]
diff_name <- diff_name[ , c(9,1,2,3,4,5,6,7,8)]

loc_diff <- merge(diff_name, name_to_loc2, by.x = "name", by.y = "name_only")
loc_diff <- loc_diff[, c(10,1:9)]

###Key next step: importing PPAg ###

setwd("per_locus_credsets/")

file_list <- list.files()

for (file in file_list){

  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }

  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }

}

setwd("~/Oxford/RealScripts/credible_sets/data")
per_locus_credset <- dataset
per_locus_credset$full_loc <- paste(per_locus_credset$Chr, per_locus_credset$Pos, sep=":")

#Merge PPA with diff
loc_PPA_diff <- merge(per_locus_credset, loc_diff, by.x = "full_loc", by.y = "gen_loc")

#Quick test
loc_count <- loc_PPA_diff[!duplicated(loc_PPA_diff$full_loc), ]
name_count <- loc_PPA_diff[!duplicated(loc_PPA_diff$name), ]

#Plot PPA vs CNN difference
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V1)
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V2)
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V3)
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V4, xlab = "PPAg", ylab = "Stage 4 CNN Difference")
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V5, xlab = "PPAg", ylab = "Stage 5 CNN Difference")
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V6, xlab = "PPAg", ylab = "Stage 6 CNN Difference")
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V7, xlab = "PPAg", ylab = "Stage 7 CNN Difference")
plot(loc_PPA_diff$PPAg, loc_PPA_diff$V8)

#Plot stage vs stage with colour gradient based on PPA
rbPal <- colorRampPalette(c('red','blue'))
loc_PPA_diff$Col <- rbPal(10)[as.numeric(cut(loc_PPA_diff$PPAg,breaks = 10))]
plot(loc_PPA_diff$V6,loc_PPA_diff$V1,pch = 20,col = loc_PPA_diff$Col)

max_test <- loc_PPA_diff$V6[loc_PPA_diff$PPAg>0.1]
#Plotting stage4(V6)PFG vs stage8(V1)BLC
library(ggplot2)
ggplot(loc_PPA_diff, aes(x=V6,
               y=V1,
               color=PPAg))+ geom_point(aes(alpha = PPAg)) +
  scale_alpha("PPAg") + scale_alpha(range = c(0.1, 1)) + scale_colour_gradient(low = "blue", high="black")

#Plotting stage5(V5)PFG vs stage8(V1)BLC
ggplot(loc_PPA_diff, aes(x=V5,
                         y=V1,
                         color=PPAg))+ geom_point(aes(alpha = PPAg)) +
  scale_alpha("PPAg") + scale_alpha(range = c(0.1, 1)) + scale_colour_gradient(low = "blue", high="black")

#Plotting stage5(V5)PFG vs stage8(V1)BLC (improved look)
ggplot(loc_PPA_diff, aes(x = V5, y = V1, color = PPAg)) +
    geom_point(aes(alpha = PPAg)) + scale_alpha("PPAg") +
    scale_colour_gradient(low = "blue", high="black") +
    annotate("point",
            loc_PPA_diff$V5[loc_PPA_diff$PPAg > 0.2],
            loc_PPA_diff$V1[loc_PPA_diff$PPAg > 0.2]) + theme_bw() +
            xlab("Stage 5: Pancreatic Endoderm (PE)") +
            ylab("Stage 8: Beta-like cells (BLC)") +
            ggtitle("Weighted scatterplot of reference and variant loci CNN predictions")

ggplot(loc_PPA_diff, aes(x = V3, y = V1, color = PPAg)) +
  geom_point(aes(alpha = PPAg)) + scale_alpha("PPAg") +
  scale_colour_gradient(low = "blue", high="black") +
  annotate("point",
           loc_PPA_diff$V3[loc_PPA_diff$PPAg > 0.2],
           loc_PPA_diff$V1[loc_PPA_diff$PPAg > 0.2]) + theme_bw() +
  xlab("Stage 7: Endocrine-like Cells (EN)") +
  ylab("Stage 8: Beta-like cells (BLC)") +
  ggtitle("Weighted scatterplot of reference and variant loci CNN predictions")

#Find the point with the highest stage 5 diff w. PPA above given value
max(loc_PPA_diff$V5[loc_PPA_diff$PPAg > 0.2])
# The location of that locus: 1:214150821, name is PROX1_rs79687284_Known_2

#Select the four lowest V5 values above PPA threshold
inc_sortedV5 <- loc_PPA_diff[order(loc_PPA_diff$V5, decreasing =  FALSE), ]
select_incsV5 <- inc_sortedV5[inc_sortedV5$PPAg > 0.2, ]
top4_neg_V5_0.2 <- select_incsV5[1:4, ]

#Look the four up in name_to_loc
#rs4976033 5:67714246 PIK3R1_rs4976033_Novel_1
#rs9521510 13:110426871 pearso_rs4771648_Novel_2
#rs11708067 3:123065778 ADCY5_rs11708067_Known_1
#rs10974438 9:4291928 GLIS3_rs10974438_Known_1

#Determine Pearson Correlation between Stage 5 and Stage 8,and Stage 7 and Stage 8
cor(loc_PPA_diff$V5, loc_PPA_diff$V1, method = "pearson")
cor(loc_PPA_diff$V3, loc_PPA_diff$V1, method = "pearson")
#Stage5 vs Stage8 = 0.004963122, Stage7 vs Stage8 = 0.9641889
