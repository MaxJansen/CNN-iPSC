###   Quick Description   ###
#After compute_credible to preprocess the tables, this script can be used to
#Determine the distribution of results 

library(dplyr)
library(fitdistrplus)
library(logspline)
library(LambertW)
library(ggplot2)
library(reshape2)
library(qvalue)
library(ape)
library(cluster) 
library(ggdendro)
library(dendextend)
library(psych)
library(gplots)

setwd("~/Oxford/RealScripts/credible_sets/data")

cred_set_results <- read.table("credible_set_largef.txt")
name_var_seq <- read.table("unique_all_name_seq.csv", header = TRUE, sep = ",")
name_to_loc <- read.table("HRC_credset.snp_ann.txt")


#Calculate the differences for each locus and each stage
#You take alternating rows and subtract the second from the first.
###   IMPORTNANT EDIT   ### 
#I switched it, now the first row (ref) is subtracted froom the second row. This way the diff_value increases if chromatin is 
#opened by a variant.
diff <- cred_set_results[seq(2,nrow(cred_set_results),2), ] - cred_set_results[seq(1,nrow(cred_set_results),2), ] 

#Nice boxplot to compare stages
#First put the stages in chronological order
diff_order <- diff 
colnames(diff_order) <- c("BLC", "DE","EN","EP","PE","PFG","PGT","iPSC")
diff_order <- diff_order[,c(8,2,7,6,5,4,3,1)]
#sapply(diff_order, sd, na.rm = TRUE)
meltDiff <- melt(diff_order)

p <- ggplot(meltDiff, aes(factor(variable), value, fill=variable)) +
  labs(title="Boxplot of Predicted Stage differences",x="Stage", y = "Predicted Difference")
p + geom_boxplot() + scale_fill_brewer(palette="RdBu") + theme_minimal()

#Now that you have the graph see if sign different
#non-parametric so kw test
kruskal.test(value ~ variable, data = meltDiff)
#Kruskal-Wallis rank sum test
#Kruskal-Wallis chi-squared = 293.49, df = 7, p-value
#< 2.2e-16 
####Significantly different variance!###

### DO NOT NORMALIZE DATA ###


#New p-values. Don't normalize the data
z_score_final <- as.data.frame(scale(diff_order))
average_diff <- colMeans(diff_order)
p_val_final <- sapply(z_score_final, function(z){pvalue2sided=2*pnorm(-abs(z))})
qfinal <- qvalue(p = p_val_final)
hist(qfinal)


summary(qfinal)
#Call:
#  qvalue(p = p_val_final)
#
#pi0:	1	
#
#Cumulative number of significant calls:#  
#  <1e-04 <0.001 <0.01 <0.025 <0.05  <0.1     <1
#p-value     8758  13796 25453  34927 46157 64742 878232
#q-value     3990   5609  8753  10788 12890 15823 878232
#local FDR   3049   4127  6033   7326  8509 10101  24284

qvalue_final <- qfinal$qvalues
lfdr <- qfinal$lfdr
hist(lfdr)
plot(qfinal)


### Merging ###

#To compare diffs to position PPA, merge dataframes
#(Optional) Tailor name_to_loc
name_to_loc2 <- name_to_loc[, -3]
colnames(name_to_loc2) <- c("name_only", "gen_loc")

#Add names and locations to diff
diff_name <- diff_order
diff_name$name <- name_var_seq$name_only[seq(1,nrow(name_var_seq),2)]
diff_name <- diff_name[ , c(9,1,2,3,4,5,6,7,8)]
loc_diff <- merge(diff_name, name_to_loc2, by.x = "name", by.y = "name_only")
loc_diff <- loc_diff[, c(10,1:9)]

### Warning! optional branch ahead ###
#To keep whole names, do this:
full_loc_diff <- merge(diff_name, name_to_loc, by.x = "name", by.y = "V1")
full_loc_diff <- full_loc_diff[, c(10,11,1:9)]

### End of optional branch

### End of first merge ###

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


### Resume optional branch ###

#Merge PPA with full_diff
colnames(full_loc_diff)[1:2] <- c("loc","nickname")
full_PPA_diff <- merge(per_locus_credset, full_loc_diff, by.x = "full_loc", by.y = "loc")
### TODO!
#Get unique ones, remove extra cols and Select PPa threshold 
full_PPA_diff <- full_PPA_diff[unique(full_PPA_diff$name),]
full_PPA_diff <- subset(full_PPA_diff, select = -c(full_loc, Chr, Pos))
full_PPA_select <- full_PPA_diff[full_PPA_diff$PPAg > 0.1, ]
###End optional branch ###



#Select only unique rows so you get 109776 rows like diff_name
loc_PPA_diff <- loc_PPA_diff[unique(loc_PPA_diff$name),]
loc_PPA_diff <- loc_PPA_diff[,c(5:14)]
#Select PPa >0.1
loc_PPA_select <- loc_PPA_diff[loc_PPA_diff$PPAg >= 0.1, ]

#Select rows with names and pred. diff. containing FDR < 0.05
fdrselect_diffname <- diff_name[apply(qvalue_final[, ], MARGIN = 1, function(x) any(x <= 0.05)), ]
qvalue_final_df <- as.data.frame(qvalue_final)

####
fdrselect_iPSC <- diff_name[qvalue_final_df$iPSC <= 0.05, ]
fdrselect_DE <- diff_name[qvalue_final_df$DE <= 0.05, ]
fdrselect_PGT <- diff_name[qvalue_final_df$PGT <= 0.05, ]
fdrselect_PFG <- diff_name[qvalue_final_df$PFG <= 0.05, ]
fdrselect_PE <- diff_name[qvalue_final_df$PE <= 0.05, ]
fdrselect_EP <- diff_name[qvalue_final_df$EP <= 0.05, ]
fdrselect_EN <- diff_name[qvalue_final_df$EN <= 0.05, ]
fdrselect_BLC <- diff_name[qvalue_final_df$BLC <= 0.05, ]
###

#Now, from the selected fdr, subset the PPa > 0.1 from loc_PPA_select. fdrselect_nn contains nicknames and variant names
fdrselect_diffname <- subset(fdrselect_diffname, name %in% loc_PPA_select$name)
fdrselect_nn <- subset(full_loc_diff, name %in% fdrselect_diffname$name)
fdrselect_nn$total_name <- paste(fdrselect_nn$name, fdrselect_nn$nickname, sep = "  ")
fdrselect_nn <- fdrselect_nn[, -c(1:3)]
rownames(fdrselect_nn) <- fdrselect_nn$total_name
fdrselect_nn <- within(fdrselect_nn, rm(total_name))

# Significant variants for each stage #
iPSC_final_select <- subset(fdrselect_iPSC, name %in% loc_PPA_select$name)
DE_final_select <- subset(fdrselect_DE, name %in% loc_PPA_select$name)
PGT_final_select <- subset(fdrselect_PGT, name %in% loc_PPA_select$name)
PFG_final_select <- subset(fdrselect_PFG, name %in% loc_PPA_select$name)
PE_final_select <- subset(fdrselect_PE, name %in% loc_PPA_select$name)
EP_final_select <- subset(fdrselect_EP, name %in% loc_PPA_select$name)
EN_final_select <- subset(fdrselect_EN, name %in% loc_PPA_select$name)
BLC_final_select <- subset(fdrselect_BLC, name %in% loc_PPA_select$name)
# End of significant variants for each stage #

#Redo heatmap with dendrogram with the following traits:
#Use fdrselect_diffname, because PPA > 0.1 and q-value < 0.05 gives 44 variants
#Add q-value asterisks
#SNP_ID + variant ID as label
fdr_nn_matrix <- data.matrix(fdrselect_nn)
heatmap.2(fdr_nn_matrix, col = redblue, dendrogram = "row", Colv = NA, key= T, cexRow = 0.6, margins = c(4,13))



#Select PROX1 variants of interest and other one and plot:
prox_selected <- full_PPA_diff[full_PPA_diff$nickname == "PROX1_rs79687284_Known_2", ]
prox_selected <- prox_selected[,-c(1:3)]
prox_selected.m <- melt(prox_selected)

p<-ggplot(prox_selected.m, aes(x=variable, y=value, group=name, linetype=name)) +
  geom_line(aes(color=name), size =1.2)+
  geom_point(aes(color=name), size =3) +
  ylim(-0.4,0.4)
p <- p + scale_color_brewer(palette="Paired")+
  labs(title="PROX1 variants throughout development",x="Stage", y = "Predicted Difference") +
  theme_minimal(base_size = 22)
p + theme(legend.position="top")

#Select ADCY5 variants of interest and other one and plot:
adcy5_selected <- full_PPA_diff[full_PPA_diff$nickname == "ADCY5_rs11708067_Known_1", ]
adcy5_selected <- adcy5_selected[adcy5_selected$PPAg > 0.1 ,]
adcy5_selected <- adcy5_selected[,-c(1:3)]
adcy5_selected.m <- melt(adcy5_selected)

p<-ggplot(adcy5_selected.m, aes(x=variable, y=value, group=name, linetype=name)) +
  geom_line(aes(color=name), size =1.2)+
  geom_point(aes(color=name), size =3) +
  ylim(-0.4,0.4)
p <- p + scale_color_brewer(palette="Paired")+
  labs(title="ADCY5 variants throughout development",x="Stage", y = "Predicted Difference") +
  theme_minimal(base_size = 22)
p + theme(legend.position="top")


#Select HNF1A variants of interest and other one and plot:
hnf1a_selected <- full_PPA_diff[full_PPA_diff$nickname == "HNF1A_rs56348580_Known_1", ]
hnf1a_selected <- hnf1a_selected[hnf1a_selected$PPAg > 0.1 ,]
hnf1a_selected <- hnf1a_selected[,-c(1:3)]
hnf1a_selected.m <- melt(hnf1a_selected)

p<-ggplot(hnf1a_selected.m, aes(x=variable, y=value, group=name, linetype=name)) +
  geom_line(aes(color=name), size =1.2)+
  geom_point(aes(color=name), size =3) +
  ylim(-0.4,0.4)
p <- p + scale_color_brewer(palette="Paired")+
  labs(title="HNF1A variants throughout development",x="Stage", y = "Predicted Difference") +
  theme_minimal(base_size = 22)
p + theme(legend.position="top")

#Select PPARG variants of interest and other one and plot:
pparg_selected <- full_PPA_diff[full_PPA_diff$nickname == "PPARG_rs17819328_Known_2", ]
pparg_selected <- pparg_selected[pparg_selected$PPAg > 0.1 ,]
pparg_selected <- pparg_selected[,-c(1:3)]
pparg_selected.m <- melt(pparg_selected)

p<-ggplot(pparg_selected.m, aes(x=variable, y=value, group=name, linetype=name)) +
  geom_line(aes(color=name), size =1.2)+
  geom_point(aes(color=name), size =3) +
  ylim(-0.4,0.4)
p <- p + scale_color_brewer(palette="Paired")+
  labs(title="PPARG variants throughout development",x="Stage", y = "Predicted Difference") +
  theme_minimal(base_size = 22)
p + theme(legend.position="top")


#####################   Not necessary   ###########################
#Hclust locus-name and variant name 
fdrselect_nn <- fdrselect_nn[,-c(1,3)]
fdrselect_nn$nickname <-  sapply(strsplit(as.character(fdrselect_nn$nickname), "\\_"), `[`, 1)
fdrselect_nn$nickname <- sub('[.]', '_', make.names(fdrselect_nn$nickname, unique=TRUE))
rownames(fdrselect_nn) <- fdrselect_nn$nickname
fdrselect_nn <- fdrselect_nn[,-1]


#Do hclust based on loci, nicknames can clarify
dist_diff <- dist(fdrselect_nn, method = "euclidean")
hc.cols <- hclust(dist((fdrselect_nn)))
plot(hc.cols, col = "#487AA1", col.main = "Black", main = "Predicted Difference Dendogram", col.lab = "#7C8071", 
     col.axis = "#F38630", lwd = 1, lty = 3, sub = "", axes = FALSE, xlab = "Locus", ylab = "Euclidian Distance", horiz = TRUE)

#Do hclust based on stages, no nicknames needed
dist_diff <- dist(fdrselect_diffname, method = "euclidean")
hc.cols <- hclust(dist(t(fdrselect_diffname)))
plot(hc.cols, col = "#487AA1", col.main = "Black", , main = "Predicted Difference Dendogram", col.lab = "#7C8071", 
     col.axis = "#F38630", lwd = 1, lty = 3, sub = "", axes = FALSE, xlab = "Stage", ylab = "Euclidian Distance")


###############################################################
###                                                         ###
###             End not necessary                           ###
###                                                         ###
###############################################################



setwd("~/Oxford/RealScripts/credible_sets/data")
save.image("~/Oxford/RealScripts/credible_sets/data/ideas_interpretation")

### Make all vs all scatterplots ###
rownames(diff_name) <- diff_name$name
diff_name <- diff_name[,-1]
plot(diff_name)

savehistory()

pairs.panels(diff_name, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE  # show density plots
)

########################################################################
###                                                                  ###
###   Select Monogenic DM Genes:                                     ###
###   GCK, HNF1A, HNF1B, HNF4A, KCNJ11, PPARG, WFS1                  ###
###                                                                  ###
########################################################################

GCK_df <- full_PPA_diff[grep("GCK", full_PPA_diff$nickname),]
HNF1A_df <- full_PPA_diff[grep("HNF1A", full_PPA_diff$nickname),]
HNF1B_df <- full_PPA_diff[grep("HNF1B", full_PPA_diff$nickname),]
HNF4A_df <- full_PPA_diff[grep("HNF4A", full_PPA_diff$nickname),]
KCNJ11_df <- full_PPA_diff[grep("KCNJ11", full_PPA_diff$nickname),]
PPARG_df <- full_PPA_diff[grep("PPARG", full_PPA_diff$nickname),]
WFS1_df <- full_PPA_diff[grep("WFS1", full_PPA_diff$nickname),]

mono_genes <- bind_rows(GCK_df, HNF1A_df,HNF1B_df, HNF4A_df, KCNJ11_df, PPARG_df, WFS1_df)

mono_genes_select <- mono_genes[mono_genes$PPAg > 0.1,]
mono_genes_select2 <- mono_genes[mono_genes$PPAg > 0.01,]


###Important additional stuff###


### Get q-values for predicted differences ###
qvalue_named <- cbind(as.data.frame(diff_name$name), qvalue_final)


###Make histogram of cred_set_results
colnames(cred_set_results) <- c("BLC", "DE","EN","EP","PE","PFG","PGT","iPSC")
cred_set_results <- cred_set_results[,c(8,2,7,6,5,4,3,1)]
cred_set_results2 <- melt(cred_set_results)
g <- ggplot(cred_set_results2,aes(x=value), fill = blue)
g <- g + geom_histogram()
g <- g + facet_wrap(~variable)
g +theme_minimal() 
