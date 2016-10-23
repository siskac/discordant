% File src/library/discordant/man/discordantRun.Rd

\name{discordantRun}
\alias{discordantRun}
\title{Run Discordant Algorithm}

\description{
   Runs discordant algorithm on two vectors of correlation coefficients.
}

\usage{
discordantRun(v1, v2, x, y = NULL, transform = TRUE, subsampling = FALSE, subSize = dim(x)[1], iter = 100, components = 3)
}

\arguments{
  \item{v1}{Vector of Pearson correlation coefficients in group 1}
  \item{v2}{Vector of Pearson correlation coefficients in group 2}
  \item{x}{m by n matrix where m are features and n samples}
  \item{y}{optional second m by n matrix, induces dual -omics analysis}
  \item{transform}{If TRUE v1 and v2 will be Fisher transformed}
  \item{subsampling}{If TRUE subsampling will be run}
  \item{subSize}{Indicates how many feature pairs to be used for subsampling. Default is the feature size in x}
  \item{iter}{Number of iterations for subsampling. Default is 100}
  \item{components}{Number of components in mixture model.}
}

\value{
  \item{discordPPVector}{Vector of differentially correlated posterior probabilities.}
  \item{discordPPMatrix}{Matrix of differentially correlated posterior probabilities where rows and columns reflect features}
  \item{classVector}{Vector of classes that have the highest posterior probability}
  \item{classMatrix}{Matrix of classes that have hte highest posterior probability where rows and columns reflect features}
  \item{probMatrix}{Matrix of posterior probabilities where rows are each molecular feature pair and columns are nine different classes}
  \item{loglik}{Final log likelihood}
}

\references{
 Siska C, Bowler R and Kechris K. The Discordant Method: A Novel Approach for Differential Correlation (2015), Bioinformatics. 32 (5): 690-696.
Lai Y, Zhang F, Nayak TK, Modarres R, Lee NH and McCaffrey TA. Concordant integrative gene set enrichment analysis of multiple large-scale two-sample expression data sets. (2014) BMC Genomics 15, S6.
Lai Y, Adam B-l, Podolsky R, She J-X. A mixture model approach to the tests of concordance and discordancd between two large-scale experiments with two sample groups. (2007) Bioinformatics 23, 1243-1250.
}

\author{
Charlotte Siska <siska.charlotte@gmail.edu>
}

\examples{

## load Data

data(TCGA_GBM_miRNA_microarray) # loads matrix called TCGA_GBM_miRNA_microarray
data(TCGA_GBM_transcript_microarray) # loads matrix called TCGA_GBM_transcript_microarray
print(colnames(TCGA_GBM_transcript_microarray)) # look at groups
groups <- c(rep(1,10), rep(2,20))

## DC analysis on only transcripts pairs

vectors <- createVectors(TCGA_GBM_transcript_microarray, groups = groups)
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcript_microarray)

## DC analysis on miRNA-transcript pairs

vectors <- createVectors(TCGA_GBM_transcript_microarray, TCGA_GBM_miRNA_microarray, groups = groups, cor.method = c("pearson"))
result <- discordantRun(vectors$v1, vectors$v2, TCGA_GBM_transcript_microarray, TCGA_GBM_miRNA_microarray)

}
\keyword{ model }
