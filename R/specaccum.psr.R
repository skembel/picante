#' Phylogenetic Species Richness Sample-Based Rarefaction Curve
#' 
#' Finds a sample-based rarefaction curve for phylogentic species richness for
#' a set of samples.
#' 
#' 
#' @param samp Community data matrix
#' @param tree A phylo tree object or a phylogenetic covariance matrix
#' @param permutations Number of permutations with method \code{method=
#' "random"}
#' @param method Species accumulation method, currently only \code{"random"}
#' is supported which adds samples in random order.
#' @param \dots Other parameters to functions
#' @return The function returns an object of class \code{"specaccum"} with
#' items:
#' 
#' \item{call}{ Function call. } \item{method}{ Accumulator method. }
#' \item{sites}{ Number of sites/samples. } \item{richness}{ The mean
#' phylogenetic species richness corresponding to number of sites/samples. }
#' \item{sd}{ The standard deviation of phylogenetic apecies accumulation
#' curve (or its standard error) estimated from permutations in \code{method =
#' "random"}. } \item{perm}{ Permutation results with \code{method = "random"}
#' and NULL in other cases. Each column in perm holds one permutation. }
#' @author Matthew Helmus \email{mrhelmus@@gmail.com} based on the
#' \code{vegan} package \link[vegan]{specaccum} function by Roeland Kindt and
#' Jari Oksanen.
#' @seealso \code{\link{psr}}, \code{\link[vegan]{specaccum}}
#' @references Gotelli N.J. & Colwell R.K. (2001) Quantifying biodiversity:
#' procedures and pitfalls in the measurement and comparison of species
#' richness. Ecology Letters, 4, 379-391\cr \cr Helmus M.R., Bland T.J.,
#' Williams C.K. & Ives A.R. (2007) Phylogenetic measures of biodiversity.
#' American Naturalist, 169, E68-E83
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' accum.sr<-specaccum(phylocom$sample, permutations = 100, method = "random")
#' plot(accum.sr, col="blue")
#' points(accum.sr$sites, accum.sr$richness, pch=19, col="blue")
#' 
#' accum.psr<-specaccum.psr(phylocom$sample, phylocom$phylo, permutations = 100, method = "random")
#' plot(accum.psr, add=TRUE, col = "red")
#' points(accum.psr$sites, accum.psr$richness, pch=19, col="red")
#' 
#' legend(5,5,legend=c("SR","PSR"),pch=c(19,19),col=c("blue","red"))
#' 
#' @export specaccum.psr
specaccum.psr <- function(samp, tree, permutations = 100, method = "random", ...) {
  # function adapted from the vegan package specaccum

  x <- as.matrix(samp)
  n <- nrow(x)
  p <- ncol(x)
  if (p == 1) {
    x <- t(x)
    n <- nrow(x)
    p <- ncol(x)
  }
  accumulator <- function(x, ind, tree) {
    n <- nrow(x)
    p <- ncol(x)
    xx <- x
    xx[1:n, 1:p] <- 0
    xx[apply(x[ind, ], 2, cumsum) > 0] <- 1
    PSV <- psv(xx, tree, compute.var = FALSE)
    PSV[, 1] * PSV[, 2]
  }
  METHODS <- c(
    "collector", "random", "exact", "rarefaction",
    "coleman"
  )
  method <- match.arg(method, METHODS)

  specaccum <- sdaccum <- sites <- perm <- NULL
  perm <- array(dim = c(n, permutations))
  for (i in 1:permutations)
  {
    r.x <- 0
    while (length(r.x) < n) {
      r.x <- accumulator(x, sample(n), tree)
    }
    perm[, i] <- r.x
  }
  sites <- 1:n
  specaccum <- apply(perm, 1, mean, na.rm = TRUE)
  sdaccum <- apply(perm, 1, sd, na.rm = TRUE)
  out <- list(call = match.call(), method = method, sites = sites, richness = specaccum, sd = sdaccum, perm = perm)
  class(out) <- "specaccum"
  out
}
