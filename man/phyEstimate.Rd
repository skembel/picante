\name{phyEstimate}
\alias{phyEstimate}
\alias{phyEstimateDisc}

\title{ Phylogenetic estimation of traits for unobserved taxa }
\description{
  Uses phylogenetic ancestral state reconstruction to estimate trait values for unobserved taxa.
}
\usage{
phyEstimate(phy, trait, method="pic", ...)
phyEstimateDisc(phy, trait, best.state=TRUE, cutoff=0.5, ...)
}

\arguments{
  \item{phy }{ phylo object }
  \item{trait }{ vector or data.frame containing trait values }
  \item{method }{ ancestral state estimation method used by \code{ace} (default="pic")}
  \item{best.state }{ estimate best-supported trait state for discrete variables? (default=TRUE)}  
  \item{cutoff }{ support cutoff required to declare a best.state }  
  \item{...}{ Additional arguments passed to \code{ace} }
}
\value{
 phyEstimate produces a data frame with columns:
  \item{est}{ Estimated trait value }
  \item{se}{ Standard error of estimated trait value }
 phyEstimateDisc produces a data frame with columns:
  \item{states 1..N}{ A column with statistical support is produced for each discrete trait state }
   \item{estimated.state}{ If best.state=TRUE, a column with the state with the highest support }
   \item{estimated.state.support}{ Statistical support for the state with the highest support }
}

\details{
These functions use phylogenetic ancestral state estimation to infer trait values for novel taxa on a phylogenetic tree, for continuous (\code{phyEstimate}) and discrete (\code{phyEstimateDisc}) traits.

The required input is a phylogenetic tree object plus a vector or data.frame containing estimated trait values for a subset of the taxa in the phylogenetic tree. Trait values for taxa that are present in the tree but not the trait data will be estimated using ancestral state estimation (Garland and Ives 2000). Briefly, for each taxon present in the tree but not the trait data, the phylogeny is rerooted at the most recent common ancestor of the novel taxon and the rest of the phylogeny, and the trait value of the novel taxon is estimated from the reconstructed trait value at the root of the rerooted phylogeny.

For \code{phyEstimateDisc}, the state with the highest support will be reported if argument \code{best.state=TRUE}. If the best-supported state's support is less than the specified \code{cutoff}, no best state is reported and a \code{NA} value will be returned.
}

\references{ 
T. Garland Jr., and A.R. Ives. 2000. Using the past to predict the present: confidence intervals for regression equations in phylogenetic comparative methods. American Naturalist 155:346364.

S.W. Kembel, M. Wu, J.A. Eisen, and J.L. Green. 2012. Incorporating 16S gene copy number information improves estimates of microbial diversity and abundance. PLoS Computational Biology 8(10):e1002743.

}
\author{ Steven Kembel <steve.kembel@gmail.com> }

\examples{
#generate random phylogeny
randtree <- rcoal(50)
#simulate trait evolution for a subset of taxa on phylogeny
randtraits <- sample(rTraitCont(randtree, sigma=10, root.value=100), 40)
#estimate trait values for "missing" taxa using PIC method
phyEstimate(randtree, randtraits, method="pic")
}
\keyword{univar}