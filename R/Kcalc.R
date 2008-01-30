`Kcalc` <-
function(x,phy) {
  
## Kcalc calculates Blomberg et al's K statistic for continuous
## traits evolving along a phylogeny. 
## See: Blomberg, S. P., T. Garland, and A. R. Ives. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57:717-745.
## For more information, please contact Simon Blomberg: 
## S.Blomberg1_at_uq.edu.au
##
## Kcalc wrapper by: David Ackerly, dackerly@berkeley.edu
## Further hacking by SB 7 August 2007 to do some error-checking and to handle
## nonultrametric trees.
## Further further hacking by SK Dec 2007 - minor bugfix
    
  if (length(x) != length(phy$tip.label)) stop(
              "Data vector and tree contain different numbers of taxa.")
  if (!setequal(names(x), phy$tip.label)) warning(
            "Taxon names in data vector and names of tip labels do not match.")

  if (!is.ultrametric(phy)) {
    mat <- vcv.phylo(phy) 
    vars <- diag(mat)
    weights <- varFixed(~vars)
  }
  else weights <- NULL
  
  x <- x[phy$tip.label]  
  matc <- vcv.phylo(phy, cor=TRUE) # correlation matrix
  ntax <- length(phy$tip.label)
  
# calculate "phylogenetic" mean via gls
  fit <- gls(x ~ 1,
             correlation=corSymm(matc[lower.tri(matc)],
             fixed=TRUE), weights=weights)
  ahat <- coef(fit)

#observed
  MSE <- fit$sigma^2
  MSE0 <- t(x-ahat) %*% (x - ahat)/(ntax-1)

#expected
  MSE0.MSE <- 1/(ntax-1) * (sum(diag(matc))-ntax/sum(solve(matc)))

  K <- MSE0/MSE / MSE0.MSE
  K
}

