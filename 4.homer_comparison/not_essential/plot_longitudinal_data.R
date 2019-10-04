#
#  PLOT LONGITUDINAL DAta
#

#   Data from TPM normalised counts
#   List of genes of interest
# IN: file with longitudinal data, list of elements to plot
# OUT: plot with longitudinal data
#
# Load necessary libraries
library(readr)  # to read big tables fast
library(ggplot2) # to plot
library(reshape2)  # to modify dataframes for ggplot2

###############################################################
#           VARIABLE
###############################################################

directory = "/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/"
filename = "31.01.2017.Differentiation_v2.gene.tpm.tsv"

# Data from input columns
donor = c("Ad2.1", "Ad3.1", "Neo1.1")  # original samples, here called donor

stage = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "ENstage6", "ENstage7") # 8 differentiation stages

# alternative names to plot:
stage_2 = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)

# rna-seq data
tpm = as.data.frame(read_tsv(
  paste(
    directory,
    filename,
    sep = ""
  )
))   # read tpm file for plotting longitudinal tpm data

##  ##  ##  ## CHANGES EVERY TIME  ##  ##  ##  ##
genes = c("TBP","NKX2-2","NKX6-1","GCG","INS","PROX1","PDX1")

# genes = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/10CPM_P12_S120_deepSplit2_signedHybrid_noOutliers/HOMER/HOMER_output/pink_known.txt",
#                    header = F)
# genes = as.character(unlist(genes))

# vector of genes to plot
output_name = "diff_taqman_genes"
output_type = "pdf"  # can be: "png", "tiff", "multiple_png", "multiple_tiff", "pdf"

# colour palette
diaPalette <-
  c("#C15858", "#6DA567", "#7883BA")  # Diabetologia palette

ncol = 4  # number of columns in multiple plots

########################################################


# QC and re-shaping of data for plotting

QC_and_reshaping = function() {
  #genes=sort(genes) # sort in alphabetical order
  plot_long = tpm[match(genes, tpm$GeneName), ]  # extracts from tpm data frame the rows that contain the genes of interest
  plot_long = na.omit(plot_long)              # remove NAs, in case there was not a match for a gene
  
  diff = setdiff(genes, tpm$GeneName)           # which ones in the genes' list are not in the table (probably have a different name)?
  
  #order the columns by stage
  
  nc = plot_long[, grepl("iPSC" , names(plot_long))]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage[-1])  {
    i = plot_long[, grepl(s , names(plot_long))]
    nc = cbind(nc, i)
  }
  
  plot_long = cbind(plot_long[c(1:2)], nc)
  rm(i, nc)
  
  gene_number = nrow(plot_long)   # how many genes to plot
  
  # melt data for ggplot2
  
  long = melt(plot_long, measure.vars = c(3:ncol(plot_long)))
  head(long)
  
  # rename stages and samples
  
  samples <- c(rep(donor, 8))
  
  long$variable = rep(stage_2, each = 3 * gene_number)                    # sample size times number of genes
  
  colnames(long)[which(names(long) == "variable")] <- "stage"
  long$Sample = rep(samples, each = gene_number)
  long$stage <- factor(long$stage, levels = stage_2)
  long$Sample = as.factor(long$Sample)
  long$GeneName = factor(long$GeneName, levels = genes)
  return(long)
  
} # this function takes parameters from the global environment
# Working directory:

if (output_type == "multiple_png" | output_type == "multiple_tiff") {
  ################ 1st part as QC script
  long = QC_and_reshaping()
  ###################### plot ########################
  
  
  p <-
    ggplot(data = long, aes(x = stage, y = value, group = Sample)) +
    ylab ("Expression (TPM)") +
    expand_limits(y = 0) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      size = 1,
      col = "#DCDCDC"
    ) +
    geom_line(aes(linetype = Sample, col = Sample), size = 1) +
    geom_point(size = 3, aes(shape = Sample, col = Sample)) +
    scale_color_manual(values = diaPalette) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(size = 2),
      axis.text.x = element_text(
        size = 16,
        face = "bold",
        angle = 45,
        vjust = 0.55
      ),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = 16, face = "bold.italic")
    ) +  # strip controls subfigure titles
    
    facet_wrap( ~ GeneName, scales = "free", ncol = ncol)
  
  if (output_type == "multiple_png") {
    png(
      paste(
        "/Users/Marta/Documents/WTCHG/DPhil/Plots/",
        output_name,
        ".png",
        sep = ""
      ),
      type = "cairo",
      antialias = "default",
      width = 7,
      height = 7,
      units = "in",
      res = 600,
      pointsize = 13
    )
    
    print(p)
    dev.off()
  }
  if (output_type == "multiple_tiff") {
    tiff(
      paste(
        "/Users/Marta/Documents/WTCHG/DPhil/Plots/",
        output_name,
        ".tiff",
        sep = ""
      ),
      type = "cairo",
      compression = "lzw",
      antialias = "default",
      width = 14,
      height = 18,
      units = "in",
      res = 600,
      pointsize = 13
    )
    
    print(p)
    dev.off()
  }
  
  
  
} else{
  if (output_type == "pdf") {
    long = QC_and_reshaping()  # reshape for plotting
    head(long)
    
    #create plots first
    plot_list = list()
    
    for (i in unique(long$GeneName))
    {
      long2 = long[long$GeneName == i, ] # subset each gene
      p <-
        ggplot(data = long2, aes(x = stage, y = value, group = Sample)) +
        ggtitle(unique(long2$GeneName)) +
        xlab ("Differentiation stages") +
        ylab ("Expression [TPM]") +
        expand_limits(y = 0) +
        geom_hline(
          yintercept = 0,
          linetype = "dashed",
          size = 1,
          col = "#DCDCDC"
        ) +
        geom_line(aes(linetype = Sample, col = Sample), size = 1) +
        scale_colour_manual(values = diaPalette) +  # diabetologia pallete
        geom_point(size = 3, aes(shape = Sample, col = Sample)) +
        #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(size = 2),
          axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold.italic"),
          legend.text = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 13, face = "bold")
        )
      
      plot_list[[i]] = p
      print(i)
    }
    # create progress bar
    library(tcltk)
    total = length(unique(long$GeneName))
    pb <- tkProgressBar(
      title = "progress bar",
      min = 0,
      max = total,
      width = 300
    )
    
    
    # I have to open pdf connection before looping so that it saves one gene on each page
    somePDFPath = paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/",
                        output_name,
                        ".pdf",
                        sep = "")
    pdf(file = somePDFPath)
    j = 1
    for (i in unique(long$GeneName)) {
      print(plot_list[[i]])
      Sys.sleep(0.1)
      setTkProgressBar(pb, j, label = paste(round(j / total * 100, 0),
                                            "% done"))
      j = j + 1
    }
    close(pb)
    dev.off()
    
  } else{
    for (genes in genes) {
      # individually gene by gene
      
      long = QC_and_reshaping()  # reshape for plotting
      head(long)
      p <-
        ggplot(data = long, aes(x = stage, y = value, group = Sample)) +
        ggtitle(unique(long$GeneName)) +
        xlab ("Differentiation stages") +
        ylab ("Expression [TPM]") +
        expand_limits(y = 0) +
        geom_hline(
          yintercept = 0,
          linetype = "dashed",
          size = 1,
          col = "#DCDCDC"
        ) +
        geom_line(aes(linetype = Sample, col = Sample), size = 1) +
        scale_colour_manual(values = diaPalette) +  # diabetologia pallete
        geom_point(size = 3, aes(shape = Sample, col = Sample)) +
        #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(size = 2),
          axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold.italic"),
          legend.text = element_text(size = 11, face = "bold"),
          legend.title = element_text(size = 13, face = "bold")
        )
      
      if (output_type == "png") {
        png(
          paste(
            "/Users/Marta/Documents/WTCHG/DPhil/Plots/",
            output_name,
            ".png",
            sep = ""
          ),
          type = "cairo",
          width = 8,
          height = 5,
          units = "in",
          res = 300,
          pointsize = 12
        )
        print(p)
        dev.off()
        
      }
      if (output_type == "tiff") {
        tiff(
          paste(
            "/Users/Marta/Documents/WTCHG/DPhil/Plots/",
            output_name,
            ".tiff",
            sep = ""
          ),
          type = "cairo",
          compression = "lzw",
          antialias = "default",
          width = 8,
          height = 5,
          units = "in",
          res = 1000,
          pointsize = 13
        )
        print(p)
        dev.off()
      }
    }
  }
  
}

# Mean + Standard error 3 replicates

library(plotrix)
library(dplyr)

long = QC_and_reshaping()  # reshape for plotting
head(long)
# calculating mean, sd and SEM
summary <- long[1:4] %>% 
  group_by(stage,GeneName,GeneID) %>% 
  summarise_all(funs(mean,sd,std.error))
summary$meanMinSE = summary$mean-summary$std.error
summary$meanPlusSE = summary$mean+summary$std.error

p <-
  ggplot(data = summary, aes(x = stage, y = mean,group = GeneName)) +
  ggtitle(unique(summary$GeneName)) +
  xlab ("Differentiation stages") +
  ylab ("Expression [TPM]") +
  expand_limits(y = 0) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    size = 1,
    col = "#DCDCDC"
  ) +
  geom_line( size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(x=stage,ymin=meanMinSE,ymax=meanPlusSE),
              alpha=0.4) +
  scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
  scale_fill_manual(values="#aaaaaa") + 
  theme_bw(base_size=12) +
  theme(
    strip.background = element_rect(fill="gray90", colour=FALSE),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour="gray90"),
    panel.margin = unit(1, "lines"),
    plot.title=element_text(vjust=1),
    axis.title.x=element_text(vjust=-0.5),
    axis.title.y=element_text(vjust=1))


cairo_pdf( paste(
  "/Users/Marta/Documents/WTCHG/DPhil/Plots/",
  output_name,
  "_average.pdf",
  sep = ""
),width=7,height=3, family="Arial")

print(p)
dev.off()

