# Training and testing CNNs
## General
Many architectures and hyperparameters have been tested
beforehand. Only the final hyperparameter settings and model(s) made it to this
directory. All scripts, except ['plot_predict_test.R'](./plot_predict_test.R),
run on Elder and often require GPU usage. Links to the necessary modules are
provided in the scripts.

The scripts fall in these three categories:

1. **Preprocessing** data.
  - ['Get dnase'](./preprocess/1.get_dnase). This is the ENCODE/Roadmap data
  used as negative training data. The first steps are identical to the tutorial in
  [Basset](https://help.github.com/en/articles/basic-writing-and-formatting-syntax#section-links).
  When you have the 'sample_beds.txt'-file (see tutorial steps), run: 'steps_after_sample_beds.sh'.
  - ['Removing overlap'](./preprocess/2.remove_overlap). Remove all peaks
  occurring in any of the 8 stages of iPSC samples from the negative training data.
  You will need a file containing the peaks occurring in iPSC samples. These are provided in the directory.
  If you want to make a similar file *de novo* see the scripts in 'Get dnase' and 'Final step'.
  - ['Shuffle and Zeros'](./preprocess/3.shuffle_and_zeros) Use only a shuffled
  subset of all available negative training sequences. Just run [`shuffle.sh`](./preprocess/3.shuffle_and_zeros/shuffle.sh)
  - ['Final step'](./preprocess/final_step) Append the negative training data
  to the islet iPSC sample data, and partition these in training, validation,
  and test sets.
2. **Training** the select CNN architecture on the preprocessed data.
3. **Testing** the model accuracy.

 [link](../some_locattion)


## On testing and random test and validation peaks
1. morestuff
2. morestuff
3. morestuff
