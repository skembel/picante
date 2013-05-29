\name{cor.table}
\alias{cor.table}
\title{ Table of correlations and P-values }
\description{
  Table of correlations with associated P-values and df, can be used with regular or independent contrast data
}
\usage{
cor.table(x, cor.method = c("pearson","spearman"),
    cor.type=c("standard","contrast")) 
}

\arguments{
  \item{x}{ Data frame of data points or contrasts at nodes }
  \item{cor.method}{ Correlation method (as \code{\link[stats]{cor}}) }
  \item{cor.type}{ Are data \code{standard} or independent \code{contrast} values?}
}

\value{
  \item{r}{Correlation values}
  \item{df}{Degrees of freedom}
  \item{P}{P-values}
}

\references{ Garland, T., Jr., P. H. Harvey, and A. R. Ives. 1992. Procedures for the analysis of comparative data using phylogenetically independent contrasts. Systematic Biology 41:18-32. }

\author{ Steven Kembel <steve.kembel@gmail.com> }
\keyword{univar}
