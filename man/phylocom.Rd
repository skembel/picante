\name{phylocom}
\alias{phylocom}
\docType{data}
\title{ Phylocom default data }
\description{
  Tree, community and trait data from the Phylocom 3.41 distribution
}
\usage{data(phylocom)}
\format{
A list with three elements
  \item{phylocom\$phylo}{Phylogenetic tree}
  \item{phylocom\$sample}{Community data}
  \item{phylocom\$traits}{Trait data}
}

\source{
  Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2007. Phylocom: software for the analysis of community phylogenetic structure and trait evolution. Version 3.41. http://www.phylodiversity.net/phylocom/.
}

\examples{
data(phylocom)
ses.mpd(phylocom$sample,cophenetic(phylocom$phylo))
}
\keyword{datasets}
