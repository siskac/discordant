
\name{createVectors}
\alias{createVectors}
\title{Create Pearson's correlation coefficient vectors based on bivariate data}

\description{
  Calculates correlation coefficients based on two groups of omics bivariate data. Currently, only two groups of samples can be specified. Used to make input for discordantRun().
}

\usage{
createVectors(x, y = NULL, groups, cor.method = c("spearman"))
}

\arguments{
  \item{x}{ExpressionSet of -omics data}
  \item{y}{optional second ExpressionSet of -omics data, induces dual -omics analysis}
  \item{groups}{n-length vector of 1s and 2s matching samples belonging to groups 1 and 2}
  \item{cor.method}{correlation method to measure association. Options are "spearman", "pearson", "bwmc" and "sparcc"}
}

\value{
  \item{v1}{List of correlation coefficients for group 1}
  \item{v2}{List of correlation coefficients for group 2}
}

\references{
Siska C, Bowler R and Kechris K. The Discordant Method: A Novel Approach for Differential Correlation. (2015) Bioinformatics. 32(5): 690-696.
Friedman J and Alm EJ. Inferring Correlation Networks from Genomic Survey Data. (2012) PLoS Computational Biology. 8:9, e1002687.
}

\details{
Creates vectors of correlation coefficents based on feature pairs within x or between x and y. The names of the vectors are the feature pairs taken from x and y.
}

\author{
Charlotte Siska <siska.charlotte@gmail.com>
}

\examples{

## load data
data("TCGA_GBM_miRNA_microarray") # loads matrix called TCGA_GBM_miRNA_microarray
data("TCGA_GBM_transcript_microarray") # loads matrix called TCGA_GBM_transcript_microarray
print(colnames(TCGA_GBM_transcript_microarray)) # look at groups

groups <- c(rep(1,10), rep(2,20))

# transcript-transcript pairs

vectors <- createVectors(TCGA_GBM_transcript_microarray, groups = groups, cor.method = c("pearson"))

# miRNA-transcript pairs

vectors <- createVectors(TCGA_GBM_transcript_microarray, TCGA_GBM_miRNA_microarray, groups = groups)

}

\keyword{datagen}
