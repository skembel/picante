\name{matrix2sample}
\alias{matrix2sample}

\title{ Convert community data matrix to Phylocom sample }
\description{
  Converts a community data matrix to a Phylocom database-format community sample
}
\usage{
matrix2sample(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ Community data matrix }
}

\value{
  Phylocom database-format community sample
}
\references{ Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.}
\author{ Steve Kembel <skembel@berkeley.edu> and Cam Webb <cwebb@oeb.harvard.edu> }
\examples{
data(phylocom)
matrix2sample(phylocom$sample)
}
\keyword{manip}
