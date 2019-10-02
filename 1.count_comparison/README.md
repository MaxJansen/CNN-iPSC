# Preliminary plots
## General remarks
Scripts in this directory have been used to check the number of ATAC-seq peaks
from the 8 stages of iPSC differentiation.
This directory is not essential for the next steps. It provides preliminary
scripts to count and visually compare the different filtering
methods available in the Kundaje Pipeline (applied to the same biological data).
These are:
- The original output of the pipeline.
- 1CPM filtering.
- And conservative filtering.

After assessing the test accuracies in the [next directory (training)](../2.train)
, a single best model was chosen. The best model was attained with training
on the 1CPM dataset.

## Order
1. Run the [counting script](./narrowpeaks_count.sh) on Elder, in the directory
containing the narrowpeaks-files.
2. Save this output for the three filtering methods in .csv-format,
with stages as rows, and methods as columns.
3. Run the [plotting script](./peaks_barplot.R) and save the plot. It should look
like [this](../GitHubdata/peaks_before_preprocessing.pdf)
