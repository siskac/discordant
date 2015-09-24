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

Discordant is a package for the analysis of molecular feature pairs derived from –omics data to determine if they correlate differently between phenotypic groups. Discordant uses a mixture model that “bins” molecular feature pairs based on their type of coexpression. More information on the algorithm can be found in Siska, et. al (submitted). Final output is summed up posterior probabilities of differential correlation bins. This package can be used to determine differential correlation within one –omics dataset or between two –omics datasets (provided that both –omics datasets were taken from the same samples). Also, the type of data can be any type of –omics, such as metabolomics, transcriptomic, proteomics, etc. as long as the data is continuous (numerical) rather than discrete (categorical, count).

The functions in the Discordant package provide a simple pipeline for moderate R users to determine differentially correlated pairs. The final output is a table of molecular feature pairs and their respective posterior probabilities. Functions have been written to allow flexibility for users in how they interpret results, which will be discussed further.

The Discordant method uses C code, which has been shown to compile on Linux. The C code is not able to compile on OSX Yosemite, however testing has not expanded to other operating systems.

### 2. Preliminaries

**Citing Discordant**

Discordant is originally derived from the Concordant algorithm written by Lai, et. al. When citing Discordant, please also include Lai, et. al in references.

Lai, Y., Adam, B. -l., Podolsky, R., and She, J.-X. (2007). A mixture model approach to the tests of concordance and discordance between two large-scale experiments with two-sample groups. Bioinformatics 23, 1243–1250.

Siska C., Bowler R.P and Kechris K. (2015). The Discordant Method: A Novel Approach to Differential Correlation. Bioinformatics. Submitted with Major Revisions.

**Installation**

Open directory that contains files discordant.R and discordant.C. First compile C code for R. Open a terminal window and run the following line

```
R CMD SHLIB discordant.c
```

This will create another file, discordant.so.

Next, open R in the same directory. To load the Discordant method, enter the following line of code.

```
source("discordant.R")
```

Now all functions should be loaded into R for use.

### 3. Quick workflow

**Brief Introduction**

Single –omics is when Discordant analysis is done within one –omics dataset. This means that all molecular features are analyzed to each other, rather than separating them by molecular type. This is mainly applicable to one –omics dataset, such as a single microarray experiment.

Dual -omics is when Discordant analysis is done with two -omics datasets. Molecular feature pairs analyzed are between the two -omics, i.e. transcript-protein, protein-metabolite, etc.

**Required Inputs**

`x`     
Bivariate m by n matrix where m are features and n are samples.

`y`  
Bivariate m by n matrix where m are features and n are samples. Optional, will induce dual -omics analysis. Samples must be matched with those in x.

`groups`  
vector which describes which group each sample belongs to using 1s and 2s

**Example run**

Load data into R.

```
load("TCGA_GBM_miRNASample.RData") # loads matrix called TCGA_GBM_miRNASample
load("TCGA_GBM_transcriptSample.RData") # loads matrix called TCGA_GBM_transcriptSample
```

Determine groups in omics data.

```
groups <- c(rep(1,10), rep(2,10))
```

*Single -omics analysis*

```
vectors <- createVectors(TCGA_GBM_transcriptSample, groups = groups)
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample)
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample)
```

*Dual -omics analysis*

```
vectors <- createVectors(TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample, groups = groups)
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
```

###4. Summary of Algorithm

...

###5. Outline of Analysis

**Create Correlation Vectors**

To run the Discordant algorithm correlation vectors respective to each group are necessary for input, which are easy to create using the function `createVectors`. Each correlation coefficient represents the linear correlation between two molecular features. The molecular features depend if a single -omics or dual -omics analysis has been performed. Correlation between molecular features in the same -omics dataset is single -omics, and correlation between molecular features in two different -omics datasets is dual -omics. Whether or not single -omics or dual -omics analysis is performed depends on whether one or two matrices are parameters for this function.

The other parameter is `groups`, which is a vector containing 1s and 2s that correspond to the location of samples in the column of the matrix for group 1 and group 2. For example, the control group is group 1 and the experimental group 2, and the location of samples corresponding to the two groups matches the locations of 1s and 2s in the group vector.

Single -omics
```
vectors <- createVectors(TCGA_GBM_transcriptSample, groups = groups)
```

Dual -omics
```
vectors <- createVectors(TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample, groups = groups)
```

createVectors has two outputs:

`v1`  
Correlation vector of molecular feature pairs corresponding to samples labeled 1 in group parameter.

`v2`  
Correlation vector of molecular feature pairs corresponding to samples labeled 2 in group parameter.


**Run Discordant Algorithm**

The Discordant Algorithm is in the function `discordantRun` which requires two correlation vectors and the original data. If the user wishes to generate their own correlation vector before inputting into the dataset, they can do so. However, the function will break if the dimenions of the datasets inserted do not match the correlation vector.

The posterior probability output of the Discordant algorithm are the summed DC posterior probabilities (those in the off-diagonal of the class matrix described in section 4). If the user wishes to observe the posterior probabilities differently, a matrix with the posterior probability of each class for each molecular feature pair is also available. 

Single -omics
```
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample)
```

Dual -omics
```
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
```

**Make Table to Summarize Results**

To ease the user in determining the posterior probability for each pair, the function `makeTable` was written. The only parameters required is the matrix of summed up discordant posterior probabilities from `discordantRun` and the data matrices.

Single -omics
```
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample)
```

Dual -omics
```
resultTable <- makeTable(result$discordPPMatrix, TCGA_GBM_transcriptSample, TCGA_GBM_miRNASample)
```
