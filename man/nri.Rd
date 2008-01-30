\name{nri}
\alias{nri}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Net Relatedness Index }
\description{
  Net Relatedness Index of community phylogenetic structure
}
\usage{
nri(samp, phylo.dist, null.model = c("taxa.labels", "sample.pool", "phylogeny.pool", "weighted.sample.pool"), runs = 99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samp}{ Community data matrix }
  \item{phylo.dist}{ Cophenetic phylogenetic distance matrix }
  \item{null.model}{ Choice of null models }
  \item{runs}{ Number of randomizations }
}

\value{
  Data frame of results
}
\references{ Webb et al. 2002, 2007 }
\author{ Steve Kembel <skembel@berkeley.edu> }
\section{Warning }{Weighted sample pool null ONLY works with presence-absence data} 
\seealso{ \code{\link{mpd}} }
