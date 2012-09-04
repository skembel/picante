\name{pglmm}
\alias{pglmm}
\alias{pglmm.sim}
\alias{pglmm.data}
\alias{pglmm.fit}
\title{Phylogenetic Generalized Linear Mixed Model}
\description{Fit the presence/absence of species across communities to environmental, trait and phylogenetic data with generalized linear mixed models.}
\usage{

pglmm.sim(tree,nsites=30,modelflag=1,figs=TRUE,second.env=TRUE,compscale = 1)
pglmm.data(modelflag=1,sim.dat=NULL,samp=NULL,tree=NULL,traits=NULL,env=NULL,Vcomp=NULL)
pglmm.fit(dat=NULL,Y=NULL,X=NULL,VV=NULL,sp.init=0.5,maxit=25,exitcountermax=50)

}

\arguments{
  \item{tree}{ Object of class phylo or a phylogenetic covariance matrix }
  \item{nsites}{ Number of sites to simulate }
  \item{modelflag}{ A number 1 - 5 indicating data set structure and corresponding to one of the models in Ives and Helmus (2011) }
  \item{figs}{ Generate figures }
  \item{compscale}{ compscale }
  \item{second.env}{ Simulate community data with two environmental covariates }
  \item{sim.dat}{ Object from \code{pglmm.sim} }
  \item{samp}{ Species (rows) by site (columns) community data matrix }
  \item{traits}{ Species X trait data matrix }
  \item{env}{ Site X environment data matrix }
  \item{Vcomp}{ Species X species matrix of pairwise repulsion }
  \item{dat}{ A \code{list} with \code{Y}, \code{X} and \code{VV} from \code{pglmm.sim} or similarly structured in the same way }
  \item{Y}{ The dependent variable, a (species times site) X 1 matrix of 0 and 1 }
  \item{X}{ The independent variables, a (species times site) X (trait + environment) matrix }
  \item{VV}{ A list of the covariance matrices one for each random effect }
  \item{sp.init}{ Initial values of the variances of the random effects (e.g., phylogenetic signal) }
  \item{maxit}{ \code{maxit} in \code{\link{optim}} used to estimate the GLMM }
  \item{exitcountermax}{ Number of iterations to estimate the fixed effect coefficients with penalized quasilikelihood and the variances of the random effects with restricted maximum likelihood }
}

\details{
 \emph{Phylogenetic Generalized Linear Mixed Models (PGLMM)} are generalized linear mixed model designed to test for phylogenetic patterns in community structure.
 Five models are implemented here and are designed to address (1) phylogenetic patterns in community structure,
 (2) phylogenetic variation in species sensitivities to environmental gradients among communities,
 (3) phylogenetic repulsion in which closely related species are less likely to co-occur,
 (4) trait-based variation in species sensitivities to environmental gradients,
 and (5) the combination traits and phylogeny to explain the variation in species occurrences among communities.
 Many other models can be designed by differently structuring the independent variable (\code{Y}), dependent variables (\code{X}), and covariance matrices (\code{VV}).
 This can be done by either editing the \code{pglmm.data} code or designing these objects by hand with custom code.
}

\value{
\code{pglmm.sim} returns a list with items:

\item{Vphylo}{ Phylogenetic covariance matrix }
\item{Vcomp}{ Repulsion matrix }
\item{Y}{ Community presence/absence matrix }
\item{X}{ Probabilities of a species being found in a community }
\item{u}{ Environmental gradient }
\item{bspp1}{ Species tolerances to environmental gradient 1 (phylogenetic signal) }
\item{bspp2}{  Species tolerances to environmental gradient 2 (no phylogenetic signal) }


\code{pglmm.data} returns a list with items:
\item{YY}{ Independent variable }
\item{VV}{ Covariance matrices for the random effects }
\item{XX}{ Dependent variables }

\code{pglmm.fit} returns a list with items:
\item{B}{ Coefficient estimates }
\item{B0}{ Initial estimates of the coefficients }
\item{s}{ The estimates of the scaling parameters (fitted variances) for the random effects (e.g., estimate of phylogenetic signal in community composition) }
\item{LL}{ Log likelihood of the final fitted PGLMM }
\item{flag}{ Did the estimation procedure converge? }
}

\references{Ives A.R. & Helmus M.R. (2011). Generalized linear mixed models for phylogenetic analyses of community structure. Ecological Monographs, 81, 511-525.}
\author{ Matthew Helmus <mrhelmus@gmail.com> and Anthony Ives <arives@wisc.edu> }
\seealso{\code{\link{phylostruct}}}
\note{
The estimation procedure can be slow especially for large data sets. Thus, start with a low value for \code{exitcountermax} to test if the data and model seem to be
correctly structured. Then increase \code{exitcountermax} to get convergence. It behooves the user to fully understand the structure of the covariance matrices
especially when designing complicated model structures. See Ives and Helmus (2011) for a discussion.
}
\examples{
\dontrun{
modelflag=1
sim.dat<-pglmm.sim(stree(16, "balanced"),nsites=30,modelflag=modelflag,second.env=TRUE,compscale=1)
str(sim.dat)

dat<-pglmm.data(modelflag=modelflag,sim.dat=sim.dat)
str(dat)

#A low number of iterations, maxit = 25 is probably good for most data sets, but exitcountermax may need to be increased depending on matrix sizes
out<-pglmm.fit(dat=dat,maxit=25,exitcountermax=30)
str(out)
out$s # The first row gives the estimate of phylogenetic signal and the second an estimate of how strongly species richness varies across communities. This later parameter is likely biologically uninformative for most research questions.
}
}

\keyword{univar}
