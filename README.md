# spatial-mc-correction-eeg

This repository contains functions for a multiple testing procedure designed for use with electroencephalographic (EEG) data.
The multiple comparisons correction is performed using a method based on detecting clusters of effects (statistically significant effects across adjacent electrodes).
The summed t-values from paired-samples t tests within each cluster (the cluster mass) is compared with a null permutation distribution of maximum cluster masses.
The null distribution is derived from random partitioning of the two groups/conditions in the dataset. 
Cluster masses that exceed a predefined percentile of the null distribution are counted as statistically significant.
Please note that this method is two-tailed in that the positive and negative difference cluster mass cutoffs are defined separately.

There is an additional option to specify a minimum cluster size which is declared statistically significant.

This method of multiple comparisons correction requires a MATLAB file specifying the channel neighbourhoods (i.e. which channels are adjacent to which other channels).

This function uses Student's paired t as the test statistic, however it can be used with other test statistics such as Yuen's t.
A function is supplied called yuend_ttest.
It is always good to test these functions on simulated and real null data which is similar to your own datasets which you will analyse.
