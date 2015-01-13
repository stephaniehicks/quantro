# quantro

### Why use quantro?

Quantile normalization is one of the most widely used multi-sample normalization tools for the analysis of noisy high-throughput data. Although it was originally developed for gene expression microarrays it is now used across many different high-throughput applications including RNAseq and ChIPseq. However, quantile normalization relies on assumptions about the data generation process that are not appropriate in some context. Unfortunately, no method exists to check for the appropriateness of these assumptions. 

For example in gene expression, we assume that observed differences between the distributions of each sample are due to only technical variation unrelated to biological variation. To normalize the samples, the distributions are forced to be the same. In general, this assumption is justified as only a minority of genes are expected to be differentially expressed between samples, but if the samples are expected to have a high percentage of global differences, it may not be appropriate to use quantile normalization as it may remove interesting global biological variation. 

The **quantro** R-package can be used to test a priori to the data analysis whether global normalization methods such as quantile normalization should be applied. Our method uses the raw unprocessed high-throughput data to test for global differences in the distributions across a set of groups. 

For help with the **quantro** R-package, there is a vignette available in the /vignettes folder.
  
# Installation

The R-package **quantro** can be [installed from the Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/quantro.html)
```s
source("http://bioconductor.org/biocLite.R")
biocLite("quantro")
```

After installation, the package can be loaded into R.
```s
library(quantro)
```

# Using quantro

The main function in the **quantro** package is `quantro()`.  The `quantro()` function needs two objects: (1) a data frame containing the samples to test for differences between their distributions with observations (rows) and samples (columns) (e.g. let's call it `mySamps`) and (2) a group level factor called `groupFactor` (let's call it `outcome`). This order of this factor variable must match the order of the columns in the `mySamps` object because it contains information about which group each sample is from.

To run the `quantro()` function, 
```
qtest <- quantro(object = mySamps, groupFactor = outcome)
qtest
```
Individual slots can be extracted using accessor methods:
```
summary(qtest)
quantroStat(qtest)
```

A permutation test is performed to assess the statistical significance of the test statistic `quantroStat` from `quantro()`. 

#### Elements in the output from `quantro()` include: 

Element | Description
--------|------------
`summary` | A list that contains (1) number of groups (`nGroups`), (2) total number of samples (`nTotSamples`) (3) number of samples in each group (`nSamplesinGroups`)
`anova` | ANOVA to test if the average medians of the distributions are different across groups
`MSbetween` | mean squared error between groups
`MSwithin` | mean squared error within groups
`quantroStat` | test statistic which is a ratio of the mean squared error between groups of distributions to the mean squared error within groups of distributions
`quantroStatPerm` | If `B` is not equal to 0, then a permutation test was performed to assess the statistical significance of `quantroStat`. These are the test statistics resulting from the permuted samples
`quantroPvalPerm` | If `B` is not equal to 0, then this is the $p$-value associated with the proportion of times the test statistics (`quantroStatPerm`) resulting from the permuted samples were larger than `quantroStat`



# Visualizing the results from the permutation test
There is a second function in the package called `quantroPlot()` which will plot the results from the permutation testing. The plot is a histogram of the  test statistics `quantroStatPerm` from the permuted samples from `quantro()` and the red line is the observed test statistic `quantroStat` from `quantro()`. 


```{r}
qtest <- quantro(object = mySamps, groupFactor = outcome)
quantroPlot(qtest)
```


Additional options in the `quantroPlot()` function include:

Element | Description
--------|-------------
xLab | the x-axis label
yLab | the y-axis label
mainLab | title of the histogram
binWidth | change the binwidth


# Bug reports
Report bugs as issues on the [GitHub repository](https://github.com/stephaniehicks/quantro)


# Contributors

* [Stephanie Hicks](https://github.com/stephaniehicks)
* [Rafael Irizarry](https://github.com/ririzarr)
