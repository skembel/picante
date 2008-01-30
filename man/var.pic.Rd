\name{var.pic}
\alias{var.pic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variance of phylogenetically independent contrasts }
\description{
  Calculate variance of standardized phylogenetically independent contrast values
}
\usage{
var.pic(x, phy, scaled = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Vector of trait data sorted in phy$tip.label order }
  \item{phy}{ phylo object }
  \item{scaled}{ TRUE or FALSE (PIC standardized or not) }
}

\value{
  Variance of phylogenetically independent contrasts across all nodes
}
\references{ Blomberg and Garland. 2002. Evolution. }
\author{ Steve Kembel <skembel@berkeley.edu> }

\seealso{ \code{\link{phylosignal}}}
