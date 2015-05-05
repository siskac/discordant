# discordant

INTRODUCTION

Discordant is a package for the analysis of molecular feature pairs derived from –omics data
to determine if they coexpress differently between phenotypic groups. Discordant uses a 
mixture model that “bins” pairs based on their type of coexpression. More information on 
the algorithm can be found in Siska, et. al (submitted). Final output is summed up posterior
probabilities of differential coexpression bins. This package can be used to determine 
differential coexpression within one –omics dataset or between two –omics datasets (provided
that both –omics datasets were taken from the same samples). Also, the type of data can be 
any type of –omics, such as metabolomics, transcriptomic, proteomic, etc. as long as the data
is continuous (numerical) rather than discrete (categorical, count).

The functions in the Discordant package provide a simple pipeline for moderate R users to 
determine differentially coexpressed pairs. The final output is a table of molecular feature
pairs and their respective posterior probabilities. Functions have been written to allow 
flexibility for users in how they interpret results, which will be discussed further.

~*~

REQUIRED INPUTS

All analyses require the following inputs:

data1/data2
An m by n matrix of expression/abundance values, where m is number of features and n is 
number of samples. Values should already be pre-processed and normalized respective to the 
type of –omics. Dataset should be separated by group, where data1 contains data for group 1 
and data2 contains data for group2. If running dual –omics, -omics datasets should be stacked
on top of each other.

featureNames
List of feature names in same order of m rows.

Dual -omics require the following inputs:

featureSize1/featureSize2
number of features in first –omics, and number of features in second –omics respectively.

featureNames1/featureNames2
names of features in first –omics, and names of features in second –omics respectively in same order of m rows in dataset1 and dataset2.

*~*

LOAD DISCORDANT INTO R

In R have working directory contain files discordant.R and discordant.c. Then load discordant
R code into R.

source("discordant.R")

Now all functions should be loaded into R for use.

*~*

DISCORDANT FUNCTIONS


fisherTrans

Purpose: Transforms Pearson’s correlation coefficients into z scores using Fisher’s method.

Arguments: 
rho		integer or numeric list of Pearson's correlation coefficients

Value:
z		integer or numeric list of transformed z scores


createVectors

Purpose: Creates vectors of correlation coefficients based on two groups of –omics bivariate data.

Arguments:
data1		1st group of bivariate normal data
data2		2nd group of bivariate normal data
multOmics	Boolean value indicating if single or multiple -omics is being analyzed
featureSize	Integer of feature size length of first -omics in data set. Value only 
		necessary if multOmics is TRUE.

Value:
v1		List of correlation coefficients for group 1
v2		List of correlation coefficients for group 2


discordantRun

Purpose: Runs discordant algorithm on two vectors of correlation coefficients.

Arguments:
data1           1st group of bivariate normal data
data2           2nd group of bivariate normal data
multOmics       Boolean value indicating if single or multiple -omics is being analyzed
featureSize     Integer of feature size length of first -omics in data set. Value only 
                necessary if multOmics is TRUE.

Value:
discordPPMatrix	Matrix of posterior probabilities where rows and columns reflect features
discordPPV	Vector of posterior probabilities
class		Vector of classes
probMatrix	Matrix of posterior probabilities where rows are each molecular feature pair 
		and columns are nine different classes
Convergence	Number of iterations for method to converge
loglik		Final log likelihood


makeTable

Purpose: Creates a table that where the first two columns are feature pairs and the third column is the posterior probability of discordance.

Arguments:
discordPPMatrix	Matrix of posterior probabilities taken from discordRun
featureNames1	List of feature names for first –omics analyzed (multiple –omics), or total 
		list of feature names (single –omics).
featuresNames2	List of feature names for second –omics analyzed (multiple –omics).

Value:
outMatrix	Matrix of posterior probabilities for all possible pairs.

*~*

EXAMPLE RUN

> library(MASS)
# for single -omics

> data1 <- mvrnorm(20,rep(0,20),diag(20))
> data2 <- mvrnorm(20,rep(0,20),diag(20))
> vectors <- createVectors(data1, data2, multOmics = FALSE)
> result <- discordRun(vectors$v1, vectors$v2, multOmics = FALSE, 20)
> resultsTable <- makeTable(result$discordPPMatrix, featureNames)

# for multiple –omics

> data1 <- mvrnorm(20,rep(0,20),diag(20))
> data2 <- mvrnorm(20,rep(0,20),diag(20))
> vectors <- createVectors(data1, data2, multOmics = TRUE, featureSize = 10)
> result <- discordRun(vectors$v1, vectors$v2, multOmics = TRUE, 10)
> resultsTable <- makeTable(result$discordPPMatrix, featureNames1, featureNames2)
