\name{unifrac}
\alias{unifrac}

\title{ Unweighted UniFrac distance between communities }
\description{ Calculates unweighted UniFrac, a phylogenetic beta diversity metric of the the unique (non-shared) fraction of total phylogenetic diversity (branch-length) between two communities. }
\usage{
unifrac(comm, tree)    
}

\arguments{
  \item{comm}{ Community data matrix }
  \item{tree}{ Object of class phylo - a rooted phylogeny}
}
\value{A dist object of the unweighted UniFrac distances between communities (the unique (non-shared) fraction of total phylogenetic diversity (branch-length) between two communities).}
\references{
Lozupone, C., Hamady, M., and Knight, R. 2006. UniFrac - an online tool for comparing microbial community diversity in a phylogenetic context. BMC Bioinformatics 7:371.
}
\author{ Steven Kembel <steve.kembel@gmail.com> }
\seealso{\code{\link{pd}}}
\note{
The supplied tree must be rooted. Single-species samples will be assigned a PD value equal to the distance from the root to the present.
}
\section{Warning }{
The UniFrac distance between samples will include the branch length connecting taxa in those samples and the root of the supplied tree. The root of the supplied tree may not be spanned by any taxa in the sample. If you want the root of your tree to correspond to the most recent ancestor of the taxa actually present in your samples, you should prune the tree before running \code{unifrac}:
\code{prunedTree <- prune.sample(sample,tree)}
}
\examples{
data(phylocom)
unifrac(phylocom$sample, phylocom$phylo)}
\keyword{univar}
