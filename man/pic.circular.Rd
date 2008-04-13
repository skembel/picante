\name{pic.circular}
\alias{pic.circular}
\title{Phylogenetic independent contrasts for circular data}
\description{
  Computes the phylogenetic independent contrasts for traits with circular distributions}
\usage{
pic.circular(x, phy, scaled = TRUE, var.contrasts = FALSE)}
\arguments{
  \item{x}{ A vector of class \code{circular, modulo = "2pi"} with length equal to \code{length(phy$edge.lable)}  }
  \item{phy}{ An object of class \code{phylo} }
  \item{scaled}{ logical, indicates whether the contrasts should be scaled with their expected variance (default to \code{TRUE}). }
  \item{var.contrasts}{ logical, indicates whether the expected variance of the contrasts should be returned (default to \code{FALSE}) }
}
\details{
  If x has names, these are matched to the dataframe, otherwise values are expected to be in the same order and an error is issued.  
}
\value{
  Vector of contrasts for trait \code{x} on phylogeny \code{phy}
}

\author{ Peter Cowan <pdc@berkeley.edu> based on the \code{ape} package \code{pic} function by Emmanuel Paradis}
\keyword{univar}

