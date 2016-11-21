\name{fishersTrans}
\alias{fishersTrans}
\title{Fisher Transformation of Pearson Correlation Coefficients to Z Scores}

\description{
   Transforms Pearsons correlation coefficients into z scores using Fishers method.
}

\usage{
fishersTrans(rho)
}

\arguments{
  \item{rho}{Integer or numeric list of Pearson's correlation coefficients}
}

\value{
Returns Fisher-transformed correlation coefficients
}

\references{
Fisher, R.A. (1915). "Frequency distribution of the values of the correlation coefficient in samples of an indefinitely large population". Biometrika (Biometrika Trust) 10 (4).
}

\examples{

## Create integer or list of Pearson's correlation coefficients.

library(MASS)
rhoV <- as.vector(cor(t(mvrnorm(10,rep(3,100),diag(100)))))

## Determine Fisher-Transformed z scores of rho
zV <- fishersTrans(rhoV)
}

\keyword{methods}
