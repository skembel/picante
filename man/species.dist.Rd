\name{species.dist}
\alias{species.dist}

\title{ Species co-occurrence distances}
\description{
  Compute interspecific distances based on patterns of species co-occurrence in communities.
}

\usage{
species.dist(x, metric = c("cij", "jaccard", "checkerboard", "doij"))
}

\arguments{
  \item{x}{ Community data matrix }
  \item{metric}{ Co-occurrence metric to use (see Details section for description) }
}

\value{
  A \code{dist} object with co-occurrences among all species pairs
}
\details{
Currently implemented co-occurrence measures (arguments to metric):
\describe{
    \item{cij}{ Schoener's index of co-occurrence }
    \item{jaccard}{ Jaccard index of co-occurrence }
    \item{checkerboard}{ Checkerboard index of co-occurrence }
    \item{doij}{ DOij index of co-occurrence }    
    }
}
\author{ Steven Kembel <skembel@uoregon.edu> }
\references{
Hardy, O.J. 2008. Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology 96:914-926.
}
\seealso{ \code{\link[vegan]{vegdist}} }
\keyword{univar}
