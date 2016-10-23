\name{splitMADOutlier}
\alias{splitMADOutlier}
\title{Outliers using left and right MAD}

\description{
  Identify features with outliers using left and right median absolute deviation (MAD).
}

\usage{
splitMADOutlier(mat, filter0 = TRUE, threshold = 2)
}

\arguments{
  \item{mat}{mxn matrix of -omics data, where rows are features and columns samples.}
  \item{filter0}{Option to filter out features if they have at least one 0 value. Default is TRUE.}
  \item{threshold}{Threshold of how many MADs outside the left or right median is used to determine features with outliers.}
}

\value{
  \item{mat.filtered}{Input matrix where features with outliers filtered out.}
  \item{index}{Index of features that have no outliers.}
}

\references{
Leys C, Klein O, Bernard P and Licata L. "Detecting Outliers: Do Not Use Standard Deviation Around the Mean, Use Absolute Deivation Around the Median." Journal of Experimental Social Psychology, 2013. 49(4), 764-766.
Magwene, PM, Willis JH, Kelly JK and Siepel A. "The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing." PLoS Computational Biology, 2011. 7(11), e1002255.
}

\examples{

## Simulate matrix of continuous -omics data.
data(TCGA_Breast_miRNASeq)

## Filter matrix based on outliers.
mat.filtered <- splitMADOutlier(TCGA_Breast_miRNASeq)$mat.filtered
}

\keyword{methods}
