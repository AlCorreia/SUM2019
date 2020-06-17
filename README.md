# Towards scalable and robust sum-product networks - SUM2019

This is the code for the paper "Towards scalable and robust sum-product networks" presented at 
the 13th international conference on Scalable Uncertainty Management (SUM 2019).

It is essentially an implementation of LearnSPN [1] in R, with a few important contributions:
1. **Robust estimates**: The model is capable of giving an estimate of the *robustness* of each classification. 
Roughly speaking, we see how much we can perturb the parameters of the SPN without changing the final output. 
The robustness value of each prediction quantifies this pertubation, with larger values for more robust predictions.

2. **Class-selective SPNs**: These are SPNs where class indicators come at the very top of the network. 
Another way to see this is that we have a mixture of SPNs, one for each class in the dataset. 
We show that *class-selective* SPNs not only outperform regular SPNs at classification tasks but also facilitate 
the computation of robustness estimates.

3. **Memory caches**: When running an SPN over a large dataset or repeated times 
(for instance, to compute robustness estimates) many nodes are evaluated at the same evidence many times. 
We show that having a simple memory cache at the nodes avoids unnecessary calculations and has considerable 
impact on the SPN's inference time.

---
If this material has been useful to you, please consider citing

Correia, Alvaro HC, and Cassio P. de Campos. "Towards scalable and robust sum-product networks." 
International Conference on Scalable Uncertainty Management. Springer, Cham, 2019.

```
@inbook{Correia2019Towards,
 address = {Germany},
 author = {Correia, Alvaro H.C. and P. de Campos, Cassio},
 booktitle = {Lecture Notes in Computer Science},
 doi = {10.1007/978-3-030-35514-2_31},
 editor = {Ben Amor, N. and Quost, B. and Theobald, M.},
 isbn = {978-3-030-35513-5},
 language = {English},
 pages = {409--422},
 publisher = {Springer},
 series = {Lecture Notes in Computer Science},
 title = {Towards scalable and robust sum-product networks},
 year = {2019}
}
```

## Usage

The experiments in the paper can be reproduced with the `sum2019.r` script.

### Accuracy comparison

The following command runs accuracy comparisons among regular SPNs, class-selective SPNs and XGBoost over all datasets.

`Rscript sum2019.r -e acc`

For each model, the script will save one csv file containing the information below per test sample. 
The files are named 
`<dataset name>_reg.csv` for regular SPN, 
`<dataset name>_sel.csv` for class-selective SPN,
`<dataset name>_xg.csv` for XGBoost.

- The true (y_true) and predicted (y_pred) classes;

- The difference in probability of the most likely and second most likely class according to the model (prob1-prob2);

- The robustness estimate for that sample (robustness);

- The log-probabilities of each class (logprob1, logprob2, ..., logprobm);

- The total number of node evaluations (count);

- The time taken to run inference over that sample (time);

- The number of node evaluations per loop, as many loops might be necessary to compute the robustness value (loop1, ...).

There are also two additionall csv files `<dataset name>_rgInfo.csv` and `<dataset name>_slInfo.csv` that contain information
about the architecture and total inference times of regular and class-selective SPNs, respectively.

### Memory caches evaluation

For running computing time comparisons for SPNs with and without memory caches over all datasets, just change the `-e` option to `mem`.

`Rscript sum2019.r -e mem`

This will save similar files as above, `<dataset name>_mem.csv` and `<dataset name>_nomem.csv` for SPNs without memory caches, respectively.

### Extra parameters

- `-d` Define a specific dataset to run on, instead of evaluating all 25 datasets (see the datasets folder).
Default is NULL which indicates that all datasets will be evaluated;

- `-f` Number of folds for the cross-validation. Default 5;

- `-r` Number of runs. Default 1;

- `g` The maximum height of the SPN. Defaul 1e6, effectively no height restriction;

- `t` The threshold (p-value) below which variables are considered independent. Default 0.01;

[1] Gens, Robert, and Domingos Pedro. 
"Learning the structure of sum-product networks." 
International conference on machine learning. 2013.
