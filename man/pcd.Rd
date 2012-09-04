\name{pcd}
\alias{pcd}
\alias{PCD}

\title{Phylogenetic Community Dissimilarity}
\description{Pairwise dissimilarity in phylogenetic community composition that is partitioned into a nonphylogenetic and a phylogenetic component.}
\usage{
    pcd(comm, tree, PSVmncd=NULL, PSVpool=NULL, reps=10^4)
}

\arguments{
  \item{comm}{ Community data matrix }
  \item{tree}{ Object of class phylo or a phylogenetic covariance matrix }
  \item{PSVmncd}{ Vector of null mean conditional \emph{phylogenetic species variability (\link[=psv]{PSV})} values }
  \item{PSVpool}{ The standard, unconditional \emph{PSV} calculated for the species pool }
  \item{reps}{ The number of random draws from the species pool used to produce \code{PSVmncd} }
}

\details{
 \emph{Phylogenetic community dissimilarity (PCD)} is the pairwise differences between communities derived by asking how much of the variance 
 among species in the values of a hypothetical nonselected trait in one community can be predicted by the known trait values of species in another community. 
 \emph{PCD} is partitioned into a nonphylogenetic component that reflects shared species between communities (\emph{PCDc}) 
 and a phylogenetic component that reflects the evolutionary relationships among nonshared species (\emph{PCDp}). In order to compare communities that vary 
 in species richness, the metric is standardized under the assumption that the species in communities are selected at random from the species pool. The   
 analyses here define the species pool as the list of all species in the set of communities in \code{comm}, but the species pool can be defined under 
 any hypothesis of community assembly either by manipulating the code or inputting a user defined \code{PSVmncd} and \code{PSVpool}.
}

\value{
The function returns a list with items:

\item{PCD}{ A square matrix of \emph{PCD} values }
\item{PCDc}{ A square matrix of \emph{PCDc} values }
\item{PCDp}{ A square matrix of \emph{PCDp} values }
\item{PSVmncd}{ A vector of null mean conditional \emph{PSV} values used to calculate \emph{PCD} }
\item{PSVpool}{ The unconditional PSV of the species pool used to calculate \emph{PCD} }
}

\references{Ives A.R. & Helmus M.R. (2010). Phylogenetic metrics of community similarity. The American Naturalist, 176, E128-E142.}
\author{ Anthony Ives <arives@wisc.edu> and Matthew Helmus <mrhelmus@gmail.com> }
\seealso{\code{\link{psv}}, \code{\link{phylosor}}, \code{\link{unifrac}}}
\note{
The sampling procedure used to standardize \emph{PCD} and produce \code{PSVmncd} and \code{PSVpool} can be slow.
}
\examples{
\dontrun{
data(phylocom)
pcd(phylocom$sample, phylocom$phylo)
}
}
\keyword{univar}
