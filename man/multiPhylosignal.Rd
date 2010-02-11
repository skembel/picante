\name{multiPhylosignal}
\alias{multiPhylosignal}
\title{ Calculates phylogenetic signal for data.frame of traits }
\description{
  Calculates phylogenetic signal for data.frame of traits. Traits may have missing values in which case the tree will be pruned prior to calculating phylogenetic signal for each trait.
}
\usage{
multiPhylosignal(x, phy, checkdata=TRUE, ...)
}
\arguments{
  \item{x}{ Data frame of trait data (traits in columns) with row names corresponding to tip.labels}
  \item{phy}{ phylo object }
  \item{checkdata}{ Check for match between trait and phylogeny taxa labels using \code{\link{match.phylo.data}}? (default=TRUE)}    
  \item{...}{ Additional arguments to \link{phylosignal} }
}
\value{
  Returns a data frame with phylogenetic signal results for each trait
}
\author{ Steven Kembel <skembel@uoregon.edu> }
\keyword{univar}
