\name{species.dist}
\alias{species.dist}

\title{ Species co-occurrence distances}
\description{
  Compute interspecific distances based on patterns of species co-occurrence in communities.
}

\usage{
species.dist(x, metric = c("cij", "jaccard", "checkerboard", "roij"))
}

\arguments{
  \item{x}{ Community data matrix }
  \item{metric}{ 
    \item{cij}{ Schoener's index of co-occurrence }
    \item{jaccard}{ Jaccard index of co-occurrence }
    \item{checkerboard}{ Checkerboard index of co-occurrence }
    \item{roij}{ ROij index of co-occurrence }    
  }
}

\value{
  A \code{dist} object with co-occurrences among all species pairs
}

\author{ Steve Kembel <skembel@berkeley.edu> }
\references{
Hardy, O. In press. Testing the spatial phylogenetic structure of local communities : 
statistical performances of different null models and test statistics on a locally neutral community. Journal of Ecology.
}
\seealso{ \code{\link[vegan]{vegdist}} }
\keyword{univar}
