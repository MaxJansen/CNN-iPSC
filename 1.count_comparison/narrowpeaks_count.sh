### Run this script on elder to count the number of peaks for
### each filtering method of the Kundaje pipeline.

#Go to the directory containing the narrowpeaks:
cd /well/mccarthy/users/martac/for_Max/narrowpeaks
#And the directory with the original peaks:

#Loop over all 3 files, unzip and count the lines to count narrowpeaks
for f in *; do echo "$f"; zcat "$f" | wc -l; done

#Save output as a .csv-file and open in R.
