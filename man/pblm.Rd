\name{pblm}
\alias{pblm}
\alias{pblmpredict}
\title{ Phylogenetic Bipartite Linear Model }
\description{
  Fits a linear model to the association strengths of a bipartite data set with or without phylogenetic correlation among the interacting species
}
\usage{
pblm(assocs, tree1=NULL, tree2=NULL, covars1=NULL, covars2=NULL, bootstrap=FALSE,
    nreps=10, maxit=10000, pstart=c(.5,.5))
pblmpredict(x, tree1.w.novel=NULL, tree2.w.novel=NULL, predict.originals=FALSE)

}

\arguments{
  \item{assocs}{ A matrix of association strengths among two sets of interacting species }
  \item{tree1}{ A phylo tree object or a phylogenetic covariance matrix for the rows of \code{assocs} }
  \item{tree2}{ A phylo tree object or a phylogenetic covariance matrix for the columns of \code{assocs}}
  \item{covars1}{ A matrix of covariates (e.g., traits) for the row species of \code{assocs} }
  \item{covars2}{ A matrix of covariates (e.g., traits) for the column species of \code{assocs} }
  \item{bootstrap}{ logical, bootstrap confidence intervals of the parameter estimates }
  \item{nreps}{ Number of bootstrap replicated data sets to estimate parameter CIs }
  \item{maxit}{ as in \code{\link{optim}} }
  \item{pstart}{ starting values of the two phylogenetic signal strength parameters passed to \code{\link{optim}} }
  \item{x}{ object of class \code{pblm} }
  \item{tree1.w.novel}{ A phylo tree object or a phylogenetic covariance matrix which corresponds to \code{tree1} of \code{x} with species to predict associations }
  \item{tree2.w.novel}{ A phylo tree object or a phylogenetic covariance matrix which corresponds to \code{tree2} of \code{x} with species to predict associations } 
  \item{predict.originals}{ if \code{TRUE} then the associations of each original species in the two phylogenies is predicted } 
}

\details{
 Fit a linear model with covariates using estimated generalized least squares to the association strengths between two sets of interacting species. 
 Associations can be either binary or continuous. If phylogenies of the two sets of interacting species are supplied, 
 two \emph{phyogenetic signal strength} parameters (\emph{d1} and \emph{d2}), one for each species set, based on an Ornstein-Uhlenbeck model of 
 evolution with stabilizing selection are estimated. Values of \emph{d=1} indicate no stabilizing selection and correspond to the Brownian motion model of 
 evolution; \emph{0<d<1} represents stabilizing selection; \emph{d=0} depicts the absence of phylogenetic correlation (i.e., a star phylogeny); and \emph{d>1} corresponds 
 to disruptive selection where phylogenetic signal is amplified. Confidence intervals for these and the other parameters can be estimated with 
 bootstrapping.
 \cr
 \cr
 The function \code{pblmpredict} predicts the associations of novel species following the methods given in appendix B of Ives and Godfray (2006). 
 }

\value{
  \item{MSE}{ total, full (each \emph{d} estimated), star (\emph{d=0}), and base (\emph{d=1}) mean squared errors }
  \item{signal.strength}{ two estimates of phylogenetic signal strength }
  \item{coefficients}{ estimated intercept and covariate coefficients with approximate 95 percent CIs for the three model types (full, star, base) }
  \item{CI.boot}{ 95 percent CIs for all parameters }
  \item{variates}{ matrix of model variates (can be used for plotting) }
  \item{residuals}{ matrix of residuals from the three models (full, star and base) }
  \item{predicted}{ predicted associations }
  \item{bootvalues}{ matrix of parameters estimated from the \code{nreps} bootstrap replicated data sets used to calculate CIs }
  \item{phylocovs}{ phylogenetic covariance matricies scaled by the estimated \code{d1} and \code{d2} }
  \item{cors.1}{ correlations among predicted and observed associations for species of \code{tree1}, \code{NA} if \code{predict.originals=FALSE} }
  \item{cors.2}{ correlations among predicted and observed associations for species of \code{tree2},  \code{NA} if \code{predict.originals=FALSE} }
  \item{pred.novels1}{ predicted associations for the novel speices of \code{tree1} }
  \item{pred.novels2}{ predicted associations for the novel speices of \code{tree2} }
}      



\note{Covariates that apply to both species sets (e.g., sampling site) should be supplied in the covariate matrix of the set with the most species.
\cr
\cr
Bootstrapping CIs is slow due to the function \code{\link{optim}} used to estimate the model parameters. See appendix A in Ives and Godfray (2006) 
for a discussion about this boostrapping procedure
\cr
\cr
If \code{pblmpredict=TRUE} the function does not first remove each species in turn when predicting the associations of the original species as 
is done in Ives and Godfray (2006).}
 
\references{Ives A.R. & Godfray H.C. (2006) Phylogenetic analysis of trophic associations. The American Naturalist, 168, E1-E14 \cr
\cr
Blomberg S.P., Garland T.J. & Ives A.R. (2003) Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution, 57, 717-745
}

\author{ Matthew Helmus \email{mrhelmus@gmail.com} }
\examples{
\dontrun{
#load example data from Ives & Godfray (2006)
data(IvesGodfray)

#net attack rate of parasitoid on host eq.4 in Ives and Godfray
A<-(-1*log(1-IvesGodfray$interactions[,-28]/t(IvesGodfray$interactions[28])))

# Make tips of the phylogenetic trees contemporaneous by extending tips
p<-dim(IvesGodfray$host)[1]
q<-dim(IvesGodfray$parasitoid)[1]
host.cov.scaled<-IvesGodfray$host
para.cov.scaled<-IvesGodfray$parasitoid
for (i in 1:p)
{
  host.cov.scaled[i,i]<-max(host.cov.scaled)
}
for (i in 1:q)
{
  para.cov.scaled[i,i]<-max(para.cov.scaled)
}

# scale covariance matrices (this reduces numerical problems caused by
# determinants going to infinity or zero)
  host.cov.scaled<-host.cov.scaled/(det(as.matrix(host.cov.scaled))^(1/p))
  para.cov.scaled<-para.cov.scaled/(det(as.matrix(para.cov.scaled))^(1/q))

pblm.A <- pblm(sqrt(A),tree1=host.cov.scaled,tree2=para.cov.scaled)
pblm.A$signal.strength    #compare to Ives and Godfray (2006) Table 1 Line 1
pblm.A$MSE
}
}
\keyword{univar}