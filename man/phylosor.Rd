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
\value{A distance object of the PhyloSor index of similarity between communities, the fraction of PD (branch-length) shared between two samples}
\references{Bryant, J.B., Lamanna, C., Morlon, H., Kerkhoff, A.J., Enquist, B.J., Green, J.L. 2008. Microbes on mountainsides: Contrasting elevational patterns of bacterial and plant diversity. Proceedings of the National Academy of Sciences 105 Supplement 1: 1505-11511}
\author{ Helene Morlon <morlon.helene@gmail.com> and Steven Kembel <skembel@uoregon.edu> }
\seealso{\code{\link{phylosor.rnd}}, \code{\link{pd}}}
\note{
The root of the supplied tree is included in calculations of PhyloSor. The supplied tree must be rooted and ultrametric. Single-species samples will be assigned a PD value equal to the distance from the root to the present.
}
\section{Warning }{
The phylosor of all samples will include the branch length connecting taxa in those samples and the root of the supplied tree. The root of the supplied tree may not be spanned by any taxa in the sample. If you want the root of your tree to correspond to the most recent ancestor of the taxa actually present in your sample, you should prune the tree before running \code{phylosor}:

\code{prunedTree <- prune.sample(sample,tree)}
}
\examples{
data(phylocom)
phylosor(phylocom$sample, phylocom$phylo)}
\keyword{univar}
