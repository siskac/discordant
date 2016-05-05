# User's Guide

###Contents

1. Introduction
2. Preliminaries
    + Citing Discordant  
    + Installation
3. Quick Start
    + Brief Introduction
    + Required Inputs
    + Example Run
4. Brief Summary of Algorithm
5. Outline of Analysis
    + Create Correlation Vectors
    + Run Discordant Algorithm
    + Make Table to Summarize Results

### 1. Introduction

*Discordant* is an R package that identifies pairs of features that correlate differently between phenotypic groups, with application to -omics datasets. *Discordant* uses a mixture model that “bins” molecular feature pairs based on their type of coexpression. More information on the algorithm can be found in Siska, et. al. The final output are posterior probabilities of differential correlation. This package can be used to determine differential correlation within one –omics dataset or between two –omics datasets (provided that both –omics datasets were taken from the same samples). Also, the type of data can be any type of –omics, such as metabolomics, transcriptomic, proteomics, etc. as long as the data are continuous (numerical) rather than discrete (categorical, count).

The functions in the *Discordant* package provide a simple pipeline for intermediate R users to determine differentially correlated pairs. The final output is a table of molecular feature pairs and their respective posterior probabilities. Functions have been written to allow flexibility for users in how they interpret results, which will be discussed further. Currently, the package only supports the comparison between two phenotypic groups (e.g., disease vs control, mutant vs wildtype).

From preliminary testing, *Discordant* compiles on OSX and Linux. We have not been able to test it on Windows. If you have tried to install this on windows, please let us know if it worked or any issues you ran into.

### 2. Preliminaries

**Citing Discordant**

*Discordant* is originally derived from the Concordant algorithm written by Lai, et. al. When citing Discordant, please also include Lai, et. al in references.

Lai, Y., Adam, B. -l., Podolsky, R., and She, J.-X. (2007). A mixture model approach to the tests of concordance and discordance between two large-scale experiments with two-sample groups. Bioinformatics 23, 1243–1250.

Siska C., Bowler R.P and Kechris K. (2015). The Discordant Method: A Novel Approach to Differential Correlation. Bioinformatics.

**Installation**

Download tarball `discordant_2.0.0.tar.gz`. In the same directory containing the tar ball, type

```
R CMD INSTALL discordant_2.0.0.tar.gz
```

Discordant has now been loaded into R. You can access discordant functions by using the `library()` function.

```
library(discordant)
```

Now all functions should be loaded into R for use.

### 3. Quick workflow

**Brief Introduction**

Single –omics refers to when the *Discordant* analysis is performed within one –omics dataset. This means that all molecular features are analyzed to each other, rather than separating them by molecular type. This is mainly applicable to one –omics dataset, such as a single microarray experiment.

Dual -omics refers to when the *Discordant* analysis is performed with two -omics datasets. Molecular feature pairs analyzed are between the two -omics, i.e. transcript-protein, protein-metabolite, etc.

**Required Inputs**

`x`     
m by n matrix where m are features and n are samples.

`y`  
m by n matrix where m are features and n are samples. Optional, will induce dual -omics analysis. Samples must be matched with those in x.

`groups`  
vector which describes which group each sample belongs to using 1s and 2s

**Example Run with Microarrays**

Load data into R.

```
data(TCGA_GBM_miRNA_microarray) # loads matrix called TCGA_GBM_miRNA_microarray
data(TCGA_GBM_transcript_microarray) # loads matrix called TCGA_GBM_microarray
```

Determine groups in omics data.

```
groups <- c(rep(1,10), rep(2,10))
```

*Single -omics analysis*

```
vectors <- createVectors(TCGA_GBM_transcript_microarray, groups = groups, cor.method = c("spearman"))
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcript_microarray)
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcript_microarray)
```

*Dual -omics analysis*

```
vectors <- createVectors(TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample, groups = groups, cor.method = c("pearson"))
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
```
**Example Run with Sequencing**

Load data into R.

```
data(TCGA_Breast_miRNASeq) # loads matrix called TCGA_Breast_miRNASeq
data(TCGA_Breast_RNASeq) # loads matrix called TCGA_Breast_RNASeq
```

Determine groups in omics data.

```
groups <- c(rep(1,15), rep(2,42))
```

*Single -omics analysis*

```
vectors <- createVectors(TCGA_Breast_RNASeq, groups = groups, cor.method = c("bwmc"))
result <- discordantRun(vectors$v1, vectors$v2, TCGA_Breast_RNASeq)
resultTable <- makeTable(result$discordPPMatrix, TCGA_Breast_RNASeq)
```

*Dual -omics analysis*

```
vectors <- createVectors(TCGA_Breast_RNASeq, TCGA_Breast_miRNASeq, groups = groups)
result <- discordantRun(vectors$v1, vectors$v2, TCGA__Breast_RNASeq, TCGA_Breast_miRNASeq)
resultTable <- makeTable(result$discordPPMatrix, TCGA_Breast_RNASeq, TCGA_Breast_miRNASeq)
```


###4. Summary of Algorithm

Using a three component mixture model and the EM algorithm, the model predicts if the correlation coefficients in phenotypic groups 1 and 2 for a molecular feature pair are dissimilar. The correlation coefficients are generated for all possible molecular feature pairs between -omics A and -omics B (Figure 1a) and are transformed in to z scores using Fisher's tranformation (Figure 1b). The three components are -, + and 0 which correspond respectively to a negative, positive or no correlation (Figure 1c). Molecular features that have correlation coefficients in *different* components are considered *differentially* correlated, as opposed when correlation coefficients are in the *same* component they are *equivalently* correlated.

![Algorithm Pipeline](siska_discordant_figure1.png)

Figure 1. Algorithm pipeline. a. Determine Pearson’s correlation coefficients for all A and B pairs. b. Fisher’s transformation c. Mixture model based on z scores d. Class matrix describing between group relationships e. EM Algorithm to estimate posterior probability of each class for each pair f. Identify features of –omics A and B that have high pp of DC.

The class matrix (Figure 1d) are the classes that represent all possible paired-correlation scenarios. These scenarios are based off the components in the mixture models. Molecular features that have correlation coefficients in *different* components are considered *differentially* correlated, as opposed to when correlation coefficients are in the *same* component they are *equivalently* correlated. This can be visualized in the class matrix, where the rows represent the components for group 1 and the columns represent the components for group 2. The classes on the diagonal represent equivalent correlation, and classes in the off-diagonal represent differential correlation.

After running the EM algorithm, we have 9 posterior probabilities for each molecular feature pair (Figure 1e) that correspond the the 9 classes in the class matrix. Since we want to summarize the probability that the molecular feature pair is differentially correlated, we sum the posterior probabilities representing the off-diagonal classes in the class matrix (Figure 1f).


###5. Outline of Analysis

**Create Correlation Vectors**

To run the *Discordant* algorithm correlation vectors respective to each group are necessary for input, which are easy to create using the function `createVectors`. Each correlation coefficient represents the linear correlation between two molecular features. The molecular features depend if a single -omics or dual -omics analysis has been performed. Correlation between molecular features in the same -omics dataset is single -omics, and correlation between molecular features in two different -omics datasets is dual -omics. Whether or not single -omics or dual -omics analysis is performed depends on whether one or two matrices are parameters for this function.

The other parameter is `groups`, which is a vector containing 1s and 2s that correspond to the location of samples in the column of the matrix for group 1 and group 2. For example, the control group is group 1 and the experimental group 2, and the location of samples corresponding to the two groups matches the locations of 1s and 2s in the group vector.

Single -omics
```
vectors <- createVectors(TCGA_GBM_transcriptSample, groups = groups)
```

Dual -omics
```
vectors <- createVectors(TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample, groups = groups)
```

**Correlation Metric**

We also have included different options for correlation metrics. This argument is called `cor.method` and its default is `"spearman"`. Other options are `"pearson"`, `"bwmc"` and `"sparcc"`. For information and comparison of Spearman, Pearson and biweight midcorrelation (bwmc) please read this paper by <a href = "http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-328">Song et al</a>.

The algorithm for SparCC was introduced by <a href = "http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687">Friedman et al</a> and is available online at <a href = "https://bitbucket.org/yonatanf/sparcc">bitbucket</a>. We use R code written by <a href = "https://github.com/huayingfang/CCLasso">Huaying Fang </a>.

createVectors has two outputs:

`v1`  
Correlation vector of molecular feature pairs corresponding to samples labeled 1 in group parameter.

`v2`  
Correlation vector of molecular feature pairs corresponding to samples labeled 2 in group parameter.


**Run *Discordant* Algorithm**

The *Discordant* Algorithm is in the function `discordantRun` requires two correlation vectors and the original data. If the user wishes to generate their own correlation vector before inputting into the dataset, they can do so. However, the function will break if the dimenions of the datasets inserted do not match the correlation vector.

The posterior probability output of the *Discordant* algorithm are the summed DC posterior probabilities (those in the off-diagonal of the class matrix described in Section 4). If the user wishes to observe the posterior probabilities differently, a matrix with the posterior probability of each class for each molecular feature pair is also available. 

Single -omics
```
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample)
```

Dual -omics
```
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
```

There are now optional arguments to make Discordant more flexible. There are two use subsampling in the EM algorithm and extend the mixture model from 3 to 5 components.

***Subsampling***

Subsampling is when independent feature pairs are drawn, ran through the EM algorithm to estimate parameters for a number of iterations, and then these paramters are used to maximize posterior probabilities for all feature pairs. There are several arguments introduced so the subsampling option can be run to the user's satisfaction. This option was introduced to make the Discordant method to run faster and also solve the independence assumption. Of course, it has its own set of issues which are explained in Siska, et al (submitted).

The argument `subsampling` must be set to `TRUE` for subsampling to be used. The number of independent feature pairs to be subsampled is determined by the argument `subSize` whose default value is the number of rows in `x`. The number of independent feature pairs must be less or equal to the number of features in `x` and `y`. The number of iterations to be run is set by the argument `iter`, whose default value is 100.

Example:

```
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample, subsampling = TRUE, iter = 200, subSize = 20)
```

***3 to 5 Components in Mixture Model***

Having 5 components instead of 3 in the mixture model allows the identification of feature pairs that have elevated differential correlation, or when there are associations in both groups in the same direction but one is more extreme. While this option introduces a new type of differential correlation, it does run longer and has less power than the 3-component mixture model.

The argument to use a 5-component mixture model instead of a 3-component model is `components` set to `5`. The default is to run the 3-component mixture model.

Example:

```
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample, components = 5)
```


**Make Table to Summarize Results**

To ease the user in determining the posterior probability for each pair, the function `makeTable` was included. The only parameters required is the matrix of summed up discordant posterior probabilities from `discordantRun` and the data matrices.

Single -omics
```
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample)
```

Dual -omics
```
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
```
