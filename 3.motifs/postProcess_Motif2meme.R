#!/usr/bin/Rscript

#This script adapts the .meme output of the MoVRs_Motif2meme.R script.
#That script converts homer-files to the .meme format, but Basset requires different tabs and spaces.
#Also, the E-value, has to be a numerical value, not NA.

evalueFormat <- function(inFile, outFile) {
  file_w_NA <- file(inFile, open = 'r')
  line <- readLines(file_w_NA)
  nuc_lines <- grepl("^[0-9]", line)
  nuc_vals <- line[nuc_lines]
  nuc_vals <- paste(" ", nuc_vals)
  nuc_vals <- gsub("*\t", "\t ", nuc_vals)
  line[nuc_lines] <- nuc_vals
  line <- paste(line, "\t", sep = "")
  line <- gsub("E= NA", "E= 0", line)
  fileConn <- file(outFile)
  writeLines(line, fileConn)
  close(fileConn)
}

#parse command line argument and run this function

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Need input and output file name, e-val = NA ===> e-val = 0\n",
       call. = TRUE)
} else {
  inFile = args[1]
  outFile = args[2]
}

evalueFormat(inFile, outFile)
