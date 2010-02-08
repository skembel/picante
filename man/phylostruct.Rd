\name{phylostruct}
\alias{phylostruct}
\title{ Permutations to Test for Phylogenetic Signal in Community Composition }
\description{ Randomize sample/community data matrices to create null distributions of given metrics
}
\usage{
phylostruct(samp, tree, env=NULL, metric=c("psv","psr","pse","psc","sppregs"), null.model=c("frequency",
            "richness","independentswap","trialswap"), runs=100, it=1000, alpha=0.05, fam="binomial")
}
\arguments{
  \item{samp}{ community data matrix, species as columns, communities as rows }
  \item{tree}{ phylo tree object or a phylogenetic covariance matrix }
  \item{env}{ environmental data matrix }
  \item{metric}{ if \code{metric="psv"}, \code{"psr"}, \code{"pse"}, or \code{"psc"} compares the observed mean of the respective metric to a null distribution at a given alpha; if \code{metric="sppregs"} compares the three correlations produced by \code{\link{sppregs}} to null distributions }
  \item{null.model}{ permutation procedure used to create the null distribution, see \code{\link{randomizeMatrix}}}
  \item{runs}{ the number of permutations to create the distribution, a rule of thumb is (number of communities)/alpha  }
  \item{it}{ the number of swaps for the independent and trial-swap null models, see \code{\link{randomizeMatrix}} }
  \item{alpha}{ probability value to compare the observed mean/correlations to a null distribution }
  \item{fam}{ as in \code{\link{sppregs}} }
}

\details{The function creates null distributions for the \code{\link{psd}} set of metrics and for the correlations of \code{\link{sppregs}} from observed community data sets.}

\value{
\item{metric}{ metric used }
\item{null.model}{ permutation used }
\item{runs}{ number of permutations }
\item{it}{ number of swaps if applicable }
\item{obs}{ observed mean value of a particular metric or the three observed correlations from \code{\link{sppregs}} }
\item{mean.null}{ mean(s) of the null distribution(s) }
\item{quantiles.null}{  quantiles of the null distribution(s) to compare to \code{obs}; determined by \code{alpha}}
\item{phylo.structure}{ if \code{obs} less than (alpha/2), \code{phylo.structure="underdispersed"}; if \code{obs} greater than (1-alpha/2), \code{phylo.structure="overdispersed"}; otherwise \code{phylo.structure="random"} and NULL if \code{metric="sppregs"}}
\item{nulls}{ null values of the distribution(s)}
}

\references{ Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007a) Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
\cr
\cr
Helmus M.R., Savage K., Diebel M.W., Maxted J.T. & Ives A.R. (2007b) Separating the determinants of phylogenetic community structure. Ecology Letters, 10, 917-925
\cr
\cr
Gotelli N.J. (2000) Null model analysis of species co-occurrence patterns. Ecology, 81, 2606-2621}

\author{ Matthew Helmus \email{mrhelmus@gmail.com} }
\seealso{ \code{\link{psd}} ,\code{\link{sppregs}}, \code{\link{randomizeMatrix}} }

\keyword{univar}