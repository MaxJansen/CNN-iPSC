### This script gets the useful information out of the html file.
### First cut out the json part. Between first <script> and next <\script>
### Get the var data = {X} part and save this as a .json file.

# Install some dependencies and load json file
install.packages("json")
library(jsonlite)
setwd("~/Oxford 2.0/Scripts/CNN_project/Data/motif_original_random")
tomtom_data <- fromJSON("tomtom.json")

#Get the necessary tables from this large data file. You will need:
# queries. These are the actual filters
# targets. These are the biological TF motifs
queries <- tomtom_data$queries
targets <- tomtom_data$targets
all_matches <- tomtom_data$all_matches