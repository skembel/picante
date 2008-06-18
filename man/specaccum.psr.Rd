\name{specaccum.psr}
\alias{specaccum.psr}

\title{ Phylogenetic Species Richness Sample-Based Rarefaction Curve }
\description{
  Finds a sample-based rarefaction curve for phylogentic species richness for a set of samples.
}
\usage{
specaccum.psr(samp, tree, permutations = 100, method = "random", ...)
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{tree}{ A phylo tree object or a phylogenetic covariance matrix}
  \item{permutations}{ Number of permutations with method \code{method= "random"} }
  \item{method}{ Species accumulation method, currently only \code{"random"} is supported which adds samples in random order. }
  \item{\dots}{ Other parameters to functions }
}

\value{
The function returns an object of class \code{"specaccum"} with items:

\item{call}{ Function call. }
\item{method}{ Accumulator method. }
\item{sites}{ Number of sites/samples. }
\item{richness}{ The mean phylogenetic species richness corresponding to number of sites/samples. }
\item{sd}{ The standard deviation of phylogenetic apecies accumulation curve (or its standard error) estimated from permutations in \code{method = "random"}. }
\item{perm}{ Permutation results with \code{method = "random"} and NULL in other cases. Each column in perm holds one permutation. }
}

\references{ Gotelli N.J. & Colwell R.K. (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology Letters, 4, 379-391\cr
              \cr Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007) Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83 }
\author{ Matthew Helmus \email{mrhelmus@gmail.com} based on the \code{vegan} package \link[vegan]{specaccum} function by Roeland Kindt and Jari Oksanen. }
\seealso{ \code{\link{psr}}, \code{\link[vegan]{specaccum}}}

\keyword{univar}
                                                                                                                 