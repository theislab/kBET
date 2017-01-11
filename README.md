# knnB

An R package to provide a test for batch effects in high-dimensional single-cell RNA sequencing data. This test is based on Pearson's chi^2 test. However, we should state a few more words on the actual test.

## Usage of the *knnB* function:

```R
#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
batch.estimate <- knnB(data, batch)
```
*knnB* creates (if `plot =TRUE`) a boxplot of the *knnB* rejection rates (for neighbourhoods and randomly chosen subsets of size *k*) and *knnB* returns a list with two parts:

* *summary*: summarizes the test results (with 95% confidence interval)
* *results*: the p-values of all tested samples 

## Variations:

The standard implementation of *knnB* performs a *k-nearest neighbour search* (if `knn = NULL`) with a pre-defined neighbourhood size *k0*, computes an optimal neighbourhood size (`heuristics = TRUE`) and finally 10% of the samples is randomly chosen to compute the test statistics itself (repeatedly by default to derive a confidence interval, `stats = 100`). For repeated runs of *knnB*, we recommend to run the k-nearest neighbour search separately:

```R
require('FNN')
# data: a matrix (rows: samples, columns: features (genes))
k0=floor(mean(table(batch)) #neighbourhood size: mean batch size 
knn <- get.knn(data, k=k0, algorithm = 'cover_tree')

#now run knnB with pre-defined nearest neighbours.
batch.estimate <- knnB(data, batch, k = k0, knn = knn)
```


