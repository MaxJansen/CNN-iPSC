#!/usr/bin/Rscript
setwd("~/Oxford 2.0/Scripts/CNN_project/motifs")
homer_file <- file("PWMs_filtered.meme", open = 'r')
line <- readLines(homer_file)
nuc_lines <- grepl("^[0-9]", line)
nuc_vals <- line[nuc_lines]
nuc_vals <- paste(" ", nuc_vals)
nuc_vals <- gsub("*\t", "\t ", nuc_vals)
line[nuc_lines] <- nuc_vals
line <- paste(line, "\t", sep="")

fileConn<-file("PWMs_final.meme")
writeLines(line, fileConn)
close(fileConn)