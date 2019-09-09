###   Quick Description   ###
#After compute_credible to preprocess the tables, this script can be used to
#Determine the distribution of results 

library(dplyr)
library(fitdistrplus)
library(logspline)
library(LambertW)
setwd("~/Oxford/RealScripts/credible_sets/data")

cred_set_results <- read.table("credible_set_sep.txt")
name_var_seq <- read.table("unique_all_name_seq.csv", header = TRUE, sep = ",")
name_to_loc <- read.table("HRC_credset.snp_ann.txt")

#Calculate the differences for each locus and each stage
#You take alternating rows and subtract the second from the first.
#The first row is the 
diff <- cred_set_results[seq(1,nrow(cred_set_results),2), ] - cred_set_results[seq(2,nrow(cred_set_results),2), ]


#Plot a histogram of the differences and fit normal distribution
m<-mean(diff$V1)
std<-sqrt(var(diff$V1))
h <- hist(diff$V1, density=50, breaks=1000, freq=TRUE, 
          xlab="Predicted difference in chromatin openness ", ylim = NULL, 
          main="Normal Curve over Histogram")
xfit <- seq(min(diff$V1), max(diff$V1), length = 40) 
yfit <- dnorm(xfit, mean = mean(diff$V1), sd = sd(diff$V1)) 
yfit <- yfit * diff(h$mids[1:2]) * length(diff$V1) 
lines(xfit, yfit, col = "black", lwd = 2)

#QQplot to be sure, and indeed... Not normally distributed.
qqnorm(cred_set_results[seq(1,nrow(cred_set_results),2), ]$V1); qqline(cred_set_results[seq(1,nrow(cred_set_results),2), ]$V1)

#Use descdist to find out which dist you are dealing with
descdist(diff$V1, discrete = FALSE)

#Example of z-scores using scale:
scale(diff$V1)
#SD is [1] 0.03102844

#Use LambertW for cool graphs:
test_norm(diff$V1)
#Kurtosis = 20. Too high!

#Estimate parameters to make normal:
mod.Lh <- MLE_LambertW(diff$V1, distname = "normal", type = "h")
summary(mod.Lh)

#See how it changed:
xx <- get_input(mod.Lh)
test_norm(xx)

#Rename rows for heatmap 
## still have to do this: add locus label ###

ord <- hclust( dist(fdrselect_diffname, method = "euclidean"), method = "ward.D" )$order

### OPTIONAL ###
#Prepare rownames, convert rsXXX names to known names
#clear_names <- merge(fdrselect_diffname, name_to_loc, by.x = "Variant", by.y = "V1")
### \OPTIONAL ###

#Make heatmap w. stage on x and genes on y
fdrselect_diffname$Variant <- (rownames(fdrselect_diffname))
fdrselect_diffname.m <- melt(fdrselect_diffname)
fdrselect_diffname.m$Variant <- factor( fdrselect_diffname.m$Variant, levels = rownames(fdrselect_diffname)[ord])
fdrselect_diffname.m$variable <- factor( fdrselect_diffname.m$variable, levels = colnames(fdrselect_diffname)[1:8] )

ggplot( fdrselect_diffname.m, aes(variable, Variant) ) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = ("blue"), high = ("red"))


########################################################################
###                                                                  ###    
### Save all tables in special directory and go on with second merge:###
###                        PPa                                       ###
###                                                                  ###
########################################################################



###Some graphs ###
#Some first picks, do manually
#HNF1A_rs56348580_Known_1, q-value 0.0959408, PPa: 0.11259

colnames(mono_genes_select[, 5:12])
mono_genes_select[mono_genes_select$name == "rs11065397", 5:12]
df_rs11065397 <- data.frame(Stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC" ), 
                            Activity = c(-0.05624551, -0.1234163, -0.0188171, 0.09596038, 0.1119376, 0.01520568, 0.02171335, 0.0122101))
df_rs11065397$Stage <- factor(df_rs11065397$Stage, levels = df_rs11065397[["Stage"]])
ggplot(data=df_rs11065397, aes(x=Stage, y=Activity, group = 1)) +
  geom_line()+
  geom_point() + 
  labs(title="Predicted difference in chromatin openness at HNF1A rs56348580") + 
  theme_minimal()

#PPARG_rs17819328_Known_2 rs4684854, q-value PGT: 0.004437023, PPa: 0.25062
mono_genes_select[mono_genes_select$name == "rs4684854", 5:12]
df_rs4684854 <- data.frame(Stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC" ), 
                           Activity = c(-0.04945779, 0.1059107, -0.1125046, -0.07721826, -0.04273896, -0.006742465, -0.009917624, -0.003167973))
df_rs4684854$Stage <- factor(df_rs4684854$Stage, levels = df_rs4684854[["Stage"]])
ggplot(data=df_rs4684854, aes(x=Stage, y=Activity, group = 1)) +
  geom_line()+
  geom_point() + 
  labs(title="Predicted difference in chromatin openness at PPARG rs17819328") + 
  theme_minimal()
