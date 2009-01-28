\name{phylosor}
\alias{phylosor}

\title{Phylogenetic index of beta-diversity PhyloSor}
\description{Fraction of branch-length shared between two communities}
\usage{
phylosor(samp, tree)    
}

\arguments{
  \item{samp}{ Community data matrix }
  \item{tree}{ Object of class phylo}
}
\value{Distance matrix reporting the PhyloSor index of similarity between communities: fraction of branch-length shared}
\references{Bryant, J.B., Lamanna, C., Morlon, H., Kerkhoff, A.J., Enquist, B.J., Green, J.L. 2008. Microbes on mountainsides: Contrasting elevational patterns of bacterial and plant diversity. Proceedings of the National Academy of Sciences 105 Supplement 1: 1505-11511}
\author{ Helene Morlon <morlon.helene@gmail.com> and Steven Kembel <skembel@uoregon.edu> }
\seealso{ \code{\link{phylosor.rnd}}}
\examples{
data(phylocom)
phylosor(phylocom$sample, phylocom$phylo)}
\keyword{univar}
