\name{sample2matrix}
\alias{sample2matrix}

\title{ Convert Phylocom sample to community data matrix }
\description{
  Convert a Phylocom database-format sample to community data matrix.
}
\usage{
sample2matrix(x)
}

\arguments{
  \item{x}{ Phylocom sample formatted data frame, a data frame with three columns:
    \itemize{
    \item{Column 1}{    Community name}
    \item{Column 2}{    Species abundance}
    \item{Column 3}{    Species name}
    }
  }
}

\references{   Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.}
\author{ Steven Kembel <skembel@uoregon.edu> and Cam Webb <cwebb@oeb.harvard.edu> }
\keyword{IO}
