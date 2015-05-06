# discordant

##Introduction

Discordant is a package for the analysis of molecular feature pairs derived from –omics data to determine if they coexpress differently between phenotypic groups. Discordant uses a mixture model that “bins” pairs based on their type of coexpression. More information on the algorithm can be found in Siska, et. al (submitted). Final output is summed up posterior
probabilities of differential coexpression bins. This package can be used to determine differential coexpression within one –omics dataset or between two –omics datasets (provided that both –omics datasets were taken from the same samples). Also, the type of data can be any type of –omics, such as metabolomics, transcriptomic, proteomic, etc. as long as the data
is continuous (numerical) rather than discrete (categorical, count).

The functions in the Discordant package provide a simple pipeline for moderate R users to determine differentially coexpressed pairs. The final output is a table of molecular feature pairs and their respective posterior probabilities. Functions have been written to allow flexibility for users in how they interpret results, which will be discussed further.

The Discordant method uses C code, which has been shown to compile on Linux. The C code is not able to compile on OSX Yosemite, however testing has not expanded to other operating systems.

##Loading Discordant into R

Open directory that contains files discordant.R and discordant.C. First compile C code for R. Open a terminal window and run the following line.

```
R CMD SHLIB discordant.c
```

This will create another file, discordant.so.

Next, open R in the same directory. To load the Discordant method, enter the following line of code.

```
source("discordant.R")
```

Now all functions should be loaded into R for use.

##Required Inputs

Required Inputs for all analyses:

Input                       | Description
----------------------------|------------
data1/data2                 | An m by n matrix of expression/abundance values, where m is number of features and n is number of samples. Values should already be pre-processed and normalized respective to the type of –omics. Dataset should be separated by group, where data1 contains data for group 1 and data2 contains data for group2. If running dual –omics, -omics datasets should be stacked on top of each other.
featureNames                | List of feature names in same order of m rows.

Required inputs for Dual -Omics:

Input                       | Description
----------------------------|------------
featureSize1/featureSize2   | Number of features in first –omics, and number of features in second –omics respectively.
featureNames1/featureNames2 | Names of features in first –omics, and names of features in second –omics respectively in same order of m rows in dataset1 and dataset2.

##Discordant Functions

###fisherTrans

Purpose: Transforms Pearson’s correlation coefficients into z scores using Fisher’s method.

Arguments                   | Description
----------------------------|------------
rho		            | Integer or numeric list of Pearson's correlation coefficients

Value                       | Description
----------------------------|------------
z		            | Integer or numeric list of transformed z scores

###createVectors

Purpose: Creates vectors of correlation coefficients based on two groups of –omics bivariate data.

Arguments                   | Description
----------------------------|------------
data1                       | 1st group of bivariate normal data
data2                       | 2nd group of bivariate normal data
multOmics	            | Boolean value indicating if single or multiple -omics is being analyzed
featureSize	            | Integer of feature size length of first -omics in data set. Value only necessary if multOmics is TRUE.

Value                       | Description
----------------------------|------------
v1                          | List of correlation coefficients for group 1
v2                          | List of correlation coefficients for group 2

###discordantRun

Purpose: Runs discordant algorithm on two vectors of correlation coefficients.

Arguments                   | Description
----------------------------|------------
data1                       | 1st group of bivariate normal data
data2                       | 2nd group of bivariate normal data
multOmics                   | Boolean value indicating if single or multiple -omics is being analyzed
featureSize                 | Integer of feature size length of first -omics in data set. Value only necessary if multOmics is TRUE.

Value                       | Description
----------------------------|------------
discordPPMatrix             | Matrix of posterior probabilities where rows and columns reflect features
class                       | Vector of classes
probMatrix                  | Matrix of posterior probabilities where rows are each molecular feature pair and columns are nine different classes
Convergence                 | Number of iterations for method to converge
loglik                      | Final log likelihood

###makeTable

Purpose: Creates a table that where the first two columns are feature pairs and the third column is the posterior probability of discordance.

Arguments                   | Description
----------------------------|------------
discordPPMatrix             | Matrix of posterior probabilities taken from discordRun
featureNames1               | List of feature names for first –omics analyzed (multiple –omics), or total list of feature names (single –omics).
featuresNames2              |List of feature names for second –omics analyzed (multiple –omics).

Value                       | Description
----------------------------|------------
outMatrix                   | Matrix of posterior probabilities for all possible pairs.

##Example Run

Simulations

```
library(MASS)

# for single -omics

data1 <- mvrnorm(20,rep(0,20),diag(20))
data2 <- mvrnorm(20,rep(0,20),diag(20))
featureNames <- 1:20

vectors <- createVectors(data1, data2, multOmics = FALSE)
result <- discordantRun(vectors$v1, vectors$v2, multOmics = FALSE, transform = TRUE, 20)
resultsTable <- makeTable(result$discordPPMatrix, multOmics = FALSE, featureNames)

# for multiple –omics

data1 <- mvrnorm(20,rep(0,20),diag(20))
data2 <- mvrnorm(20,rep(0,20),diag(20))
featureNames1 <- 1:10
featureNames2 <- 11:20

vectors <- createVectors(data1, data2, multOmics = TRUE, featureSize = 10)
result <- discordantRun(vectors$v1, vectors$v2, multOmics = TRUE, transform = TRUE, 10)
resultsTable <- makeTable(result$discordPPMatrix, multOmics = TRUE, featureNames1, featureNames2)
```

Real Data

```
load("TCGA_GBM_miRNASample.RData") # loads matrix called miRNASampleMatrix
load("TCGA_GBM_transcriptSample.RData") # loads matrix called transSampleMatrix
print(colnames(transSampleMatrix)) # look at groups
group1 <- 1:10
group2 <- 11:20

# DC analysis on only transcripts pairs

featureNames <- rownames(transSampleMatrix)

vectors <- createVectors(transSampleMatrix[,group1], transSampleMatrix[,group2], multOmics = FALSE)
result <- discordantRun(vectors$v1, vectors$v2, multOmics = FALSE, transform = TRUE, 20)
resultsTable <- makeTable(result$discordPPMatrix, multOmics = FALSE, featureNames)


# DC analysis on miRNA-transcript pairs

featureNames1 <- rownames(microSampleMatrix)
featureNames2 <- rownames(transSampleMatrix)
data <- rbind(microSampleMatrix, transSampleMatrix)

vectors <- createVectors(data[,group1], data[,group2], multOmics = TRUE, featureSize = dim(microSampleMatrix)[1])
result <- discordantRun(vectors$v1, vectors$v2, multOmics = TRUE, transform = TRUE, dim(microSampleMatrix)[1])
resultsTable <- makeTable(result$discordPPMatrix, multOmics = TRUE, featureNames1, featureNames2)
