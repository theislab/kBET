---
title: "kBET short introduction"
author: "Maren Büttner"
date: "9/18/2017"
---

# kBET - k-nearest neighbour batch effect test

The R package provides a test for batch effects in high-dimensional single-cell RNA sequencing data. It evaluates the accordance of replicates based on Pearson's $\chi^2$ test. First, the algorithm creates k-nearest neighbour matrix and choses 10% of the samples to check the batch label distribution in its neighbourhood. If the local batch label distribution is sufficiently similar to the global batch label distribution, the $\chi^2$-test does not reject the null hypothesis (that is "all batches are well-mixed"). The neighbourhood size k is fixed for all tests. Next, the test returns a binary result for each of the tested samples. Finally, the result of kBET is the average test rejection rate. The lower the test result, the less bias is introduced by the batch effect. kBET is very sensitive to any kind of bias. If kBET returns an average rejection rate of 1 for your batch-corrected data, you may also consider to compute the average silhouette width and PCA-based batch-effect measures to explore the degree of the batch effect. 
Learn more about kBET and batch effect correction in our [publication](https://www.nature.com/articles/s41592-018-0254-1).

## Installation

Installation should take less than 5 min. 

### Via Github and devtools

If you want to install the package directly from Github, I recommend to use the `devtools` package.

```R
library(devtools)
install_github('theislab/kBET')
```

### Manually

Please download the package as `zip` archive and install it via

```R
install.packages('kBET.zip', repos = NULL, type = 'source')
```

## Usage of the *kBET* function:

```R
#data: a matrix (rows: cells or other observations, columns: features (genes); will be transposed if necessary)
#batch: vector or factor with batch label of each cell/observation; length has to match the size of the corresponding data dimension  
batch.estimate <- kBET(data, batch)
```
*kBET* creates (if `plot = TRUE`) a boxplot of the *kBET* rejection rates (for neighbourhoods and randomly chosen subsets of size *k*) and *kBET* returns a list with several parts:

* *summary*: summarizes the test results (with 95% confidence interval)
* *results*: the p-values of all tested samples 
* *average.pval*: an average over all p-values of the tested samples 
* *stats*: the results for each of `n_repeat` runs - they can be used to reproduce the boxplot that is returned by kBET
* *params*: the parameters used in kBET
* *outsider*: samples without mutual nearest neighbour, their batch labels and a p-value whether their batch label composition varies from the global batch label frequencies

For a single-cell RNAseq dataset with less than 1,000 samples, the estimated run time is less than 2 minutes. 

### Plot *kBET*'s rejection rate

By default (`plot = TRUE`), *kBET* returns a boxplot of observed and expected rejection rates for a data set. You might want to turn off the display of these plots and create them elsewhere. *kBET* returns all information that is needed in the `stats` part of the results. 

``` R
library(ggplot2)
batch.estimate <- kBET(data, batch, plot=FALSE)
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))
g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
     labs(x='Test', y='Rejection rate',title='kBET test results') +
     theme_bw() +  
     scale_y_continuous(limits=c(0,1))
```

## Variations:

The standard implementation of *kBET* performs a *k-nearest neighbour search* (if `knn = NULL`) with a pre-defined neighbourhood size *k0*, computes an optimal neighbourhood size (`heuristics = TRUE`) and finally 10% of the samples is randomly chosen to compute the test statistics itself (repeatedly by default to derive a confidence interval, `n_repeat = 100`). For repeated runs of *kBET*, we recommend to run the k-nearest neighbour search separately:

```R
require('FNN')
# data: a matrix (rows: samples, columns: features (genes))
k0=floor(mean(table(batch))) #neighbourhood size: mean batch size 
knn <- get.knn(data, k=k0, algorithm = 'cover_tree')

#now run kBET with pre-defined nearest neighbours.
batch.estimate <- kBET(data, batch, k = k0, knn = knn)
```
It must be noted that the `get.knn` function from the `FNN` package initializes a variable with `n * k` entries, where `n` is the sample size and `k` is the neighbourhood size. If `n * k > 2^31`, the `get.knn` aborts the k-nearest neighbour search. The initial neighbourhood size in kBET (`k0`) is ~ `1/4* mean(batch size)`, which can be already too large for example for mass cytometry data. In such cases, we recommend to subsample the data. 

## Subsampling:

Currently (July 2019), *kBET* operates only on dense matrices, which results in memory issues for large datasets. Furthermore, k-nearest neighbour search with `FNN` is limited (see above). We recommend to subsample in these cases. We have thought of several options. One option is to subsample the data irrespective of the substructure:

```
#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
subset_size <- 0.1 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
batch.estimate <- kBET(data[subset_id,], batch[subset_id])
```

In case of differently sized batches, one should consider *stratified sampling* in order to keep more samples from smaller batches. 

The second option of subsampling is to take into account the substructure of the data (*i.e.* clusters). We observed that the batch label frequencies may vary in clusters. For example, such changes are due to inter-individual variability, or due to targeted population enrichment in some batches (e.g. by FACS), in contrast to unbiased cell sampling. In these cases, we compute the rejection rates for each cluster separately and average the results afterwards. 

```
#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
#clusters: vector or factor with cluster label of each cell 

kBET_result_list <- list()
sum_kBET <- 0
for (cluster_level in unique(clusters)){
   batch_tmp <- batch[clusters == cluster_level]
   data_tmp <- data[clusters == cluster_level,]
   kBET_tmp <- kBET(df=data_tmp, batch=batch_tmp, plot=FALSE)
   kBET_result_list[[cluster_level]] <- kBET_tmp
   sum_kBET <- sum_kBET + kBET_tmp$summary$kBET.observed[1]
}

#averaging
mean_kBET = sum_kBET/length(unique(clusters))
```

## Compute a silhouette width and PCA-based measure:

```R
#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
pca.data <- prcomp(data, center=TRUE) #compute PCA representation of the data
batch.silhouette <- batch_sil(pca.data, batch)
batch.pca <- pcRegression(pca.data, batch)
```
For a single-cell RNAseq dataset with less than 1,000 samples, the estimated run time is less than 2 minutes. 

