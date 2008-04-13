\name{cor.table}
\alias{cor.table}
\title{ Table of correlations and P-values }
\description{
  Table of correlations with associated P-values and df, can be used with regular or independent contrast data
}
\usage{
cor.table(x, cor.method = "pearson", cor.type=c("standard","contrast")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Data frame of data points or contrasts at nodes }
  \item{cor.method}{ Correlation method (as \code{\link{cor}}) }
  \item{cor.type}{ Are data \code{standard} or independent \code{contrast} values?}
}

\value{
  \item{r}{Correlation values}
  \item{df}{Degrees of freedom}
  \item{P}{P-values}
}

\references{ Garland et al. 1992 }

\author{ Steve Kembel <skembel@berkeley.edu> }
\keyword{univar}
