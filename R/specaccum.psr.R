specaccum.psr<-function (samp, tree, permutations = 100, method = "random", ...)
{

#function adapted from the vegan package specaccum

  x <- as.matrix(samp)
  n <- nrow(x)
  p <- ncol(x)
  if (p == 1)
  {
    x <- t(x)
    n <- nrow(x)
    p <- ncol(x)
  }
  accumulator <- function(x,ind,tree)
  {
    n <- nrow(x)
    p <- ncol(x)
    xx<-x
    xx[1:n,1:p]<-0
    xx[apply(x[ind, ], 2, cumsum)>0]<-1
    PSV<-psv(xx,tree,compute.var=FALSE)
    PSV[,1]*PSV[,2]
  }
  METHODS <- c("collector", "random", "exact", "rarefaction",
        "coleman")
    method <- match.arg(method, METHODS)

  specaccum <- sdaccum <- sites <- perm <- NULL
  perm <- array(dim = c(n, permutations))
  for (i in 1:permutations)
  {
    r.x=0
    while(length(r.x)<n){r.x <- accumulator(x, sample(n),tree)}
    perm[, i]<-r.x
  }
  sites <- 1:n
  specaccum <- apply(perm, 1, mean, na.rm=TRUE)
  sdaccum <- apply(perm, 1, sd, na.rm=TRUE)
  out <- list(call = match.call(), method = method, sites = sites, richness = specaccum, sd = sdaccum, perm = perm)
  class(out) <- "specaccum"
  out
}