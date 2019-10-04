# Preprocessing, training and testing CNNs
## General
Many architectures and hyperparameters have been tested
beforehand. Only the final hyperparameter settings and model(s) made it to this
directory. All scripts, except ['plot_predict_test.R'](./plot_predict_test.R),
run on Elder and often require GPU usage. Links to the necessary modules and packages are
provided in the scripts.

## Order
The scripts fall in these three categories:

1. **Preprocessing** data.
  - ['Get dnase'](./preprocess/1.get_dnase). This is the ENCODE/Roadmap data
  used as negative training data. The first steps are identical to the tutorial in
  [Basset](https://github.com/davek44/Basset/blob/master/tutorials/prepare_compendium.ipynb).
  When you have the 'sample_beds.txt'-file (see tutorial steps), run: 'steps_after_sample_beds.sh'
  to get 1000 bp peaks of these negative set sequences.
  - ['Removing overlap'](./preprocess/2.remove_overlap). Remove all peaks
  occurring in any of the 8 stages of iPSC samples from the negative training data.
  You will need a file containing the peaks occurring in iPSC samples. These are provided in the directory.
  If you want to make a similar file *de novo* see the scripts in 'Get dnase' and 'Final step'.
  - ['Shuffle and Zeros'](./preprocess/3.shuffle_and_zeros) Use only a shuffled
  subset of all available negative training sequences. Just run [`shuffle.sh`](./preprocess/3.shuffle_and_zeros/shuffle.sh)
  - ['Final step'](./preprocess/final_step) Append the negative training data
  to the islet iPSC sample data, and partition these in training, validation,
  and test sets. Also, clips off the negative dataset column, because the algorithm
  shouldn't optimise predicting these.
2. **Training** the select CNN architecture on the preprocessed data.
 - Once you have the `1CPM_islets.h5` file, you can use that as training data.
 - The model [hyperparameter settings](./train/filt21_params.txt) are in: `/well/mccarthy/users/maxlouis/oxford2/CNN_project/better_train_CNN/filt21_params.txt`
 on Elder.
 - The training (and testing) happens in: `/well/mccarthy/users/maxlouis/oxford2/CNN_project/better_train_CNN/1CPM_random`.
 - This is quite straightforward. The necessary paths for training data and settings
 are in the `Basset_train.sh` script. You can run it on the GPU by running
 `Basset_submit_train.sh` script. Just make sure the GPUtype is identical in both scripts
(you can switch from `p100` to `k80`, depending on which one is free).
3. **Testing** the model accuracy.
For testing, do not use the `Basset_test.py` or the `Basset_test.sh` script. It is the
standard testing script in Basset, but it does not give you the
freedom to plot all 8 stage lines in one ROC or PR plot.
Instead, use `Basset_predict.py`.

4. **Plotting** the test results (locally). Use the `iter1.test.txt` output as the test set
predictions. `plot_predict_test.R` takes these predictions and compares them to
`final_test_set_act.txt` to calculate the ROC and PR and plot the respective curves.
