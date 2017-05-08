# kBET - k-nearest neighbour batch effect test

An R package to provide a test for batch effects in high-dimensional single-cell RNA sequencing data. It evaluates the accordance of replicates based on Pearson's chi^2 test. First, the algorithm creates k-nearest neighbour matrix and choses 10% of the samples to check the batch label distribution in its neighbourhood. If the local batch label distribution is sufficiently similar to the global batch label distribution, the chi^2-test does not reject the null hypothesis (that is "all batches are well-mixed"). The neighbourhood size k is fixed for all tests. Next, the test returns a binary result for each of the tested samples. Finally, the result of kBET is the average test rejection rate. The lower the test result, the less bias is introduced by the batch effect.  

## Usage of the *kBET* function:

```R
#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
batch.estimate <- kBET(data, batch)
```
*kBET* creates (if `plot =TRUE`) a boxplot of the *kBET* rejection rates (for neighbourhoods and randomly chosen subsets of size *k*) and *kBET* returns a list with two parts:

* *summary*: summarizes the test results (with 95% confidence interval)
* *results*: the p-values of all tested samples 

## Variations:

The standard implementation of *kBET* performs a *k-nearest neighbour search* (if `knn = NULL`) with a pre-defined neighbourhood size *k0*, computes an optimal neighbourhood size (`heuristics = TRUE`) and finally 10% of the samples is randomly chosen to compute the test statistics itself (repeatedly by default to derive a confidence interval, `stats = 100`). For repeated runs of *kBET*, we recommend to run the k-nearest neighbour search separately:

```R
require('FNN')
# data: a matrix (rows: samples, columns: features (genes))
k0=floor(mean(table(batch))) #neighbourhood size: mean batch size 
knn <- get.knn(data, k=k0, algorithm = 'cover_tree')

#now run kBET with pre-defined nearest neighbours.
batch.estimate <- kBET(data, batch, k = k0, knn = knn)
```


