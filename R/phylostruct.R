##' Permutations to Test for Phylogenetic Signal in Community Composition
##' 
##' Randomize sample/community data matrices to create null distributions of
##' given metrics
##' 
##' The function creates null distributions for the \code{\link{psd}} set of
##' metrics and for the correlations of \code{\link{sppregs}} from observed
##' community data sets.
##' 
##' @param samp community data matrix, species as columns, communities as rows
##' @param tree phylo tree object or a phylogenetic covariance matrix
##' @param env environmental data matrix
##' @param metric if \code{metric="psv"}, \code{"psr"}, \code{"pse"}, or
##' \code{"psc"} compares the observed mean of the respective metric to a null
##' distribution at a given alpha; if \code{metric="sppregs"} compares the
##' three correlations produced by \code{\link{sppregs}} to null distributions
##' @param null.model permutation procedure used to create the null
##' distribution, see \code{\link{randomizeMatrix}}
##' @param runs the number of permutations to create the distribution, a rule
##' of thumb is (number of communities)/alpha
##' @param it the number of swaps for the independent and trial-swap null
##' models, see \code{\link{randomizeMatrix}}
##' @param alpha probability value to compare the observed mean/correlations to
##' a null distribution
##' @param fam as in \code{\link{sppregs}}
##' @return \item{metric}{ metric used } \item{null.model}{ permutation used }
##' \item{runs}{ number of permutations } \item{it}{ number of swaps if
##' applicable } \item{obs}{ observed mean value of a particular metric or the
##' three observed correlations from \code{\link{sppregs}} } \item{mean.null}{
##' mean(s) of the null distribution(s) } \item{quantiles.null}{ quantiles of
##' the null distribution(s) to compare to \code{obs}; determined by
##' \code{alpha}} \item{phylo.structure}{ if \code{obs} less than (alpha/2),
##' \code{phylo.structure="underdispersed"}; if \code{obs} greater than
##' (1-alpha/2), \code{phylo.structure="overdispersed"}; otherwise
##' \code{phylo.structure="random"} and NULL if \code{metric="sppregs"}}
##' \item{nulls}{ null values of the distribution(s)}
##' @author Matthew Helmus \email{mrhelmus@@gmail.com}
##' @seealso \code{\link{psd}} ,\code{\link{sppregs}},
##' \code{\link{randomizeMatrix}}
##' @references Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007a)
##' Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
##' \cr \cr Helmus M.R., Savage K., Diebel M.W., Maxted J.T. & Ives A.R.
##' (2007b) Separating the determinants of phylogenetic community structure.
##' Ecology Letters, 10, 917-925 \cr \cr Gotelli N.J. (2000) Null model
##' analysis of species co-occurrence patterns. Ecology, 81, 2606-2621
##' @keywords univar
##' @export phylostruct
phylostruct <- function(samp, tree, env = NULL, metric = c("psv", "psr", "pse", "psc", "sppregs"), null.model = c("frequency", "richness", "independentswap", "trialswap"), runs = 100, it = 1000, alpha = 0.05, fam = "binomial") {
  metric <- match.arg(metric)
  null.model <- match.arg(null.model)
  if (metric == "sppregs") {
    nulls <- t(replicate(runs, sppregs(randomizeMatrix(samp, null.model = null.model, iterations = it), env, tree, fam = fam)$correlations))
    obs <- sppregs(samp, env, tree, fam = fam)$correlations
    mean.null <- apply(nulls, 2, mean)
    quantiles.null <- t(apply(nulls, 2, quantile, probs = c(alpha / 2, 1 - (alpha / 2))))
    if ((null.model != "independentswap") && (null.model != "trialswap")) {
      it <- NA
    }
    return(list(
      metric = metric, null.model = null.model, runs = runs, it = it, obs = obs, mean.null = mean.null,
      quantiles.null = quantiles.null, phylo.structure = NULL, nulls = nulls
    ))
  } else {
    nulls <- switch(metric,
      psv = replicate(runs, mean(psv(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree, compute.var = FALSE)[, 1], na.rm = TRUE)),
      psr = replicate(runs, mean(psr(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree, compute.var = FALSE)[, 1], na.rm = TRUE)),
      pse = replicate(runs, mean(pse(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree)[, 1], na.rm = TRUE)),
      psc = replicate(runs, mean(psc(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree)[, 1], na.rm = TRUE))
    )
    quantiles.null <- quantile(nulls, probs = c(alpha / 2, 1 - (alpha / 2)))
    mean.null <- mean(nulls)
    mean.obs <- switch(metric,
      psv = mean(psv(samp, tree, compute.var = FALSE)[, 1], na.rm = TRUE),
      psr = mean(psr(samp, tree, compute.var = FALSE)[, 1], na.rm = TRUE),
      pse = mean(pse(samp, tree)[, 1], na.rm = TRUE),
      psc = mean(psc(samp, tree)[, 1], na.rm = TRUE)
    )

    if (mean.obs <= quantiles.null[1]) {
      phylo.structure <- "underdispersed"
    } else {
      if (mean.obs >= quantiles.null[2]) {
        phylo.structure <- "overdispersed"
      } else {
        phylo.structure <- "random"
      }
    }
    if ((null.model != "independentswap") && (null.model != "trialswap")) {
      it <- NA
    }
    return(list(
      metric = metric, null.model = null.model, runs = runs, it = it, mean.obs = mean.obs, mean.null = mean.null,
      quantiles.null = quantiles.null, phylo.structure = phylo.structure, null.means = nulls
    ))
  }
}
