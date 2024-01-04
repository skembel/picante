#' Correlations between species co-occurrence and phylogenetic distances
#' 
#' Calculates measures of community phylogenetic structure (correlation
#' between co-occurrence and phylogenetic distance) to patterns expected under
#' various null models
#' 
#' Currently implemented null models (arguments to null.model): \describe{
#' \item{sample.taxa.labels}{Shuffle phylogeny tip labels (only within set of
#' taxa present in community data)} \item{pool.taxa.labels}{Shuffle phylogeny
#' tip labels (across all taxa included in phylogenetic tree)}
#' \item{frequency}{Randomize community data matrix abundances within species
#' (maintains species occurence frequency)} \item{richness}{Randomize
#' community data matrix abundances within samples (maintains sample species
#' richness)} \item{independentswap}{Randomize community data matrix
#' maintaining species occurrence frequency and site richnessing using
#' independent swap} \item{trialswap}{Randomize community data matrix
#' maintaining species occurrence frequency and site richnessing using trial
#' swap} }
#' 
#' @param samp Community data matrix
#' @param phylo Phylogenetic tree
#' @param metric Metric of co-occurrence to use (see
#' \code{\link{species.dist}})
#' @param null.model Null model to use (see Details section for description)
#' @param runs Number of runs (randomizations)
#' @param ...  Additional arguments to \link{randomizeMatrix}
#' @return A list with elements: \item{obs.corr }{ Observed
#' co-occurrence/phylogenetic distance correlation} \item{obs.corr.p}{ P-value
#' of observed correlation (standard P-value for correlation coefficient, not
#' based on comparison with randomizations)} \item{obs.rank}{ Rank of observed
#' correlation vs. random} \item{runs}{ Number of runs (randomizations) }
#' \item{obs.rand.p}{ P-value of observed correlation vs. randomizations (=
#' obs.rank / (runs + 1))} \item{random.corrs}{A vector of random correlation
#' calculated for each run}
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{randomizeMatrix}}
#' @references Cavender-Bares J., D.A. Ackerly, D. Baum and F.A. Bazzaz. 2004.
#' Phylogenetic overdispersion in Floridian oak communities, American
#' Naturalist, 163(6):823-843.
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' comm.phylo.cor(phylocom$sample, phylocom$phylo, metric="cij",null.model="sample.taxa.labels")
#' @export comm.phylo.cor
`comm.phylo.cor` <-
  function(samp, phylo, metric = c("cij", "checkerboard", "jaccard", "doij"),
           null.model = c(
             "sample.taxa.labels", "pool.taxa.labels",
             "frequency", "richness", "independentswap", "trialswap"
           ),
           runs = 999, ...) {
    metric <- match.arg(metric)
    null.model <- match.arg(null.model)
    results <- list(
      "obs.corr" = NA, "obs.corr.p" = NA, "obs.rank" = NA, "runs" = runs,
      "obs.rand.p" = NA, "random.corrs" = vector(length = runs)
    )
    phylo.dist <- as.dist(cophenetic(prune.sample(samp, phylo)))
    pool.phylo.dist <- as.dist(cophenetic(phylo))
    taxa.names <- rownames(as.matrix(phylo.dist))
    samp.dist <- as.dist(as.matrix(species.dist(samp, metric))[taxa.names, taxa.names])
    results$obs.corr <- cor(phylo.dist, samp.dist, use = "pairwise")
    results$obs.corr.p <- cor.test(phylo.dist, samp.dist)$p.value
    if (null.model == "sample.taxa.labels") {
      for (run in 1:runs)
      {
        phylo.dist <- as.dist(taxaShuffle(as.matrix(phylo.dist))[taxa.names, taxa.names])
        results$random.corrs[run] <- cor(phylo.dist, samp.dist, use = "pairwise")
      }
    } else if (null.model == "pool.taxa.labels") {
      for (run in 1:runs)
      {
        phylo.dist <- as.dist(taxaShuffle(as.matrix(pool.phylo.dist))[taxa.names, taxa.names])
        results$random.corrs[run] <- cor(phylo.dist, samp.dist, use = "pairwise")
      }
    } else {
      for (run in 1:runs)
      {
        samp.dist <- species.dist(randomizeMatrix(samp, null.model, ...), metric)
        results$random.corrs[run] <- cor(phylo.dist, samp.dist, use = "pairwise")
      }
    }
    results$obs.rank <- rank(as.vector(c(results$obs.corr, results$random.corrs)))[1]
    results$obs.rand.p <- results$obs.rank / (runs + 1)
    results
  }


#' Quantile regression slopes between species co-occurrence and phylogenetic
#' distances
#' 
#' Calculates measures of community phylogenetic structure (quantile
#' regression between co-occurrence and phylogenetic distance) to patterns
#' expected under various null models
#' 
#' This function fits a quantile regression of co-occurrence versus
#' phylogenetic distances separating species, and compares observed patterns
#' to the patterns expected under some null model. The quantile regressions
#' are fit using the \code{\link[quantreg]{rq}} function from the
#' \code{quantreg} package.
#' 
#' Currently implemented null models (arguments to null.model): \describe{
#' \item{sample.taxa.labels}{Shuffle phylogeny tip labels (only within set of
#' taxa present in community data)} \item{pool.taxa.labels}{Shuffle phylogeny
#' tip labels (across all taxa included in phylogenetic tree)}
#' \item{frequency}{Randomize community data matrix abundances within species
#' (maintains species occurence frequency)} \item{richness}{Randomize
#' community data matrix abundances within samples (maintains sample species
#' richness)} \item{independentswap}{Randomize community data matrix
#' maintaining species occurrence frequency and site richnessing using
#' independent swap} \item{trialswap}{Randomize community data matrix
#' maintaining species occurrence frequency and site richnessing using trial
#' swap} }
#' 
#' @param samp Community data matrix
#' @param phylo Phylogenetic tree
#' @param metric Metric of co-occurrence to use (see
#' \code{\link{species.dist}})
#' @param null.model Null model to use (see Details section for description)
#' @param quant Quantile of slope to be fit (using \code{\link[quantreg]{rq}})
#' @param runs Number of runs (randomizations)
#' @param show.plot Option to display a plot of co-occurrence versus
#' phylogenetic distance with quantile regression slope fit
#' @param ...  Additional arguments to \link{randomizeMatrix}
#' @return A list with elements: \item{obs.qr.intercept }{ Observed
#' co-occurrence/phylogenetic distance quantile regression intercept}
#' \item{obs.qr.slope }{ Observed co-occurrence/phylogenetic distance quantile
#' regression slope} \item{obs.qr.slope.p}{ P-value of observed quantile
#' regression slope significance versus null model (calculated based on
#' comparison with randomizations)} \item{obs.rank}{ Rank of observed quantile
#' regression slope vs. random} \item{runs}{ Number of runs (randomizations) }
#' \item{random.qr.slopes}{A vector of quantile regression slopes calculated
#' for each randomization}
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{randomizeMatrix}}
#' @references Cavender-Bares J., D.A. Ackerly, D. Baum and F.A. Bazzaz. 2004.
#' Phylogenetic overdispersion in Floridian oak communities, American
#' Naturalist, 163(6):823-843. Slingsby, J. A. and G. A. Verboom. 2006.
#' Phylogenetic relatedness limits coexistence at fine spatial scales:
#' evidence from the schoenoid sedges (Cyperaceae: Schoeneae) of the Cape
#' Floristic Region, South Africa. The American Naturalist 168:14-27.
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' comm.phylo.qr(phylocom$sample, phylocom$phylo, metric="cij",
#'   null.model="sample.taxa.labels", runs=99)
#' @export comm.phylo.qr
`comm.phylo.qr` <-
  function(samp, phylo, metric = c("cij", "checkerboard", "jaccard", "doij"),
           null.model = c(
             "sample.taxa.labels", "pool.taxa.labels",
             "frequency", "richness", "independentswap", "trialswap"
           ),
           quant = 0.75, runs = 999, show.plot = FALSE, ...) {
    if (!requireNamespace("quantreg")) {
      stop("The 'quantreg' package is required to use this function.")
    }

    metric <- match.arg(metric)
    null.model <- match.arg(null.model)
    results <- list(
      "obs.qr.intercept" = NA, "obs.qr.slope" = NA, "obs.qr.slope.p" = NA, "obs.rank" = NA, "runs" = runs,
      "random.qr.slopes" = vector(length = runs)
    )
    phylo.dist <- as.dist(cophenetic(prune.sample(samp, phylo)))
    pool.phylo.dist <- as.dist(cophenetic(phylo))
    taxa.names <- rownames(as.matrix(phylo.dist))
    samp.dist <- as.dist(as.matrix(species.dist(samp, metric))[taxa.names, taxa.names])
    results$quantile <- quant
    qrres <- coef(quantreg::rq(samp.dist ~ phylo.dist, tau = quant, na.action = na.omit))
    names(qrres) <- NULL
    results$obs.qr.intercept <- qrres[1]
    results$obs.qr.slope <- qrres[2]

    if (show.plot) {
      plot(samp.dist ~ phylo.dist, xlab = "Phylogenetic distance", ylab = "Co-occurrence")
    }

    if (null.model == "sample.taxa.labels") {
      for (run in 1:runs)
      {
        phylo.dist <- as.dist(taxaShuffle(as.matrix(phylo.dist))[taxa.names, taxa.names])
        results$random.qr.slopes[run] <- coef(quantreg::rq(samp.dist ~ phylo.dist,
                                                           tau = quant,
                                                           na.action = na.omit
        ))[2]
      }
    } else if (null.model == "pool.taxa.labels") {
      for (run in 1:runs)
      {
        phylo.dist <- as.dist(taxaShuffle(as.matrix(pool.phylo.dist))[taxa.names, taxa.names])
        results$random.qr.slopes[run] <- coef(quantreg::rq(samp.dist ~ phylo.dist,
                                                           tau = quant,
                                                           na.action = na.omit
        ))[2]
      }
    } else {
      for (run in 1:runs)
      {
        samp.dist <- species.dist(randomizeMatrix(samp, null.model, ...), metric)
        results$random.qr.slopes[run] <- coef(quantreg::rq(samp.dist ~ phylo.dist,
                                                           tau = quant,
                                                           na.action = na.omit
        ))[2]
      }
    }
    results$obs.rank <- rank(as.vector(c(results$obs.qr.slope, results$random.qr.slopes)))[1]
    results$obs.qr.slope.p <- results$obs.rank / (runs + 1)

    if (show.plot) {
      abline(results$obs.qr.intercept, results$obs.qr.slope)
      legend("topleft", paste0("q", as.character(quant)))
    }

    results
  }


`taxaShuffle` <-
  function(x) {
    # TODO replace with vegan's permuted.index?
    if (!is.matrix(x)) x <- as.matrix(x)
    rand.names <- sample(rownames(x))
    rownames(x) <- rand.names
    colnames(x) <- rand.names
    return(x)
  }

`tipShuffle` <-
  function(phy) {
    phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
    return(phy)
  }


#' Mean pairwise distance
#' 
#' Calculates mean pairwise distance separating taxa in a community
#' @param samp Community data matrix
#' @param dis Interspecific distance matrix
#' @param abundance.weighted Should mean pairwise distances be weighted by
#' species abundance? (default = FALSE)
#' @return Vector of MPD values for each community
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{ses.mpd}}
#' @references Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and Systematics
#' 33:475-505.
#' @keywords univar
#' @examples
#' data(phylocom)
#' mpd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#' 
#' @export mpd
mpd <- function(samp, dis, abundance.weighted = FALSE)
{
  N <- dim(samp)[1] ## le nombre de row
  mpd <- numeric(N) ## N en liste de N objets
  for (i in 1:N) {
    sppInSample <- names(samp[i, samp[i,] > 0]) ## prend les rows avec des espèces présentes
    if (length(sppInSample) > 1) { ## si plus de 1
      sample.dis <- dis[sppInSample, sppInSample] ##les distances en tableau croisé
      if (abundance.weighted) {
        sample.weights <- t(as.matrix(samp[i, sppInSample, drop = FALSE])) %*% as.matrix(samp[i, sppInSample, drop = FALSE])
        mpd[i] <- weighted.mean(sample.dis, sample.weights)
      } else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
      }
    } else {
      mpd[i] <- NA
    }
  }
  mpd
}



#' Mean nearest taxon distance
#' 
#' Calculates MNTD (mean nearest taxon distance) for taxa in a community
#' 
#' This metric has also been referred to as MNND (mean nearest neighbour
#' distance), and the function was named \code{mnnd} in picante versions <
#' 0.7.
#' 
#' @aliases mntd mnnd
#' @param samp Community data matrix
#' @param dis Interspecific distance matrix
#' @param abundance.weighted Should mean nearest taxon distances for each
#' species be weighted by species abundance? (default = FALSE)
#' @return Vector of MNTD values for each community.
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{ses.mntd}}
#' @references Webb, C., D. Ackerly, M. McPeek, and M. Donoghue. 2002.
#' Phylogenies and community ecology. Annual Review of Ecology and Systematics
#' 33:475-505.
#' @keywords univar
#' @examples
#' data(phylocom)
#' mntd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
#' 
#' @export mntd
mntd <- function(samp, dis, abundance.weighted = FALSE) {
  N <- dim(samp)[1]
  mntd <- numeric(N)
  for (i in 1:N) {
    sppInSample <- names(samp[i, samp[i,] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample, sppInSample]
      ##
      diag(sample.dis) <- NA
      if (abundance.weighted) {
        mntds <- apply(sample.dis, 2, min, na.rm = TRUE)
        sample.weights <- samp[i, sppInSample]
        mntd[i] <- weighted.mean(mntds, sample.weights)
      } else {
        mntd[i] <- mean(apply(sample.dis, 2, min, na.rm = TRUE))
      }
    } else {
      mntd[i] <- NA
    }
  }
  mntd
}


#' Standardized effect size of MPD
#' 
#' Standardized effect size of mean pairwise distances in communities. When
#' used with a phylogenetic distance matrix, equivalent to -1 times the
#' Nearest Relative Index (NRI).
#' 
#' Currently implemented null models (arguments to null.model): \describe{
#' \item{taxa.labels}{ Shuffle distance matrix labels (across all taxa
#' included in distance matrix)} \item{richness}{ Randomize community data
#' matrix abundances within samples (maintains sample species richness)}
#' \item{frequency}{ Randomize community data matrix abundances within species
#' (maintains species occurence frequency)} \item{sample.pool}{ Randomize
#' community data matrix by drawing species from pool of species occurring in
#' at least one community (sample pool) with equal probability}
#' \item{phylogeny.pool}{ Randomize community data matrix by drawing species
#' from pool of species occurring in the distance matrix (phylogeny pool) with
#' equal probability} \item{independentswap}{ Randomize community data matrix
#' with the independent swap algorithm (Gotelli 2000) maintaining species
#' occurrence frequency and sample species richness } \item{trialswap}{
#' Randomize community data matrix with the trial-swap algorithm (Miklos &
#' Podani 2004) maintaining species occurrence frequency and sample species
#' richness } }
#' 
#' @param samp Community data matrix
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each
#' species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for
#' independent swap and trial null models)
#' @return A data frame of results for each community \item{ntaxa}{Number of
#' taxa in community} \item{mpd.obs}{Observed mpd in community}
#' \item{mpd.rand.mean}{Mean mpd in null communities}
#' \item{mpd.rand.sd}{Standard deviation of mpd in null communities}
#' \item{mpd.obs.rank}{Rank of observed mpd vs. null communities}
#' \item{mpd.obs.z}{Standardized effect size of mpd vs. null communities (=
#' (mpd.obs - mpd.rand.mean) / mpd.rand.sd, equivalent to -NRI)}
#' \item{mpd.obs.p}{P-value (quantile) of observed mpd vs. null communities (=
#' mpd.obs.rank / runs + 1)} \item{runs}{Number of randomizations}
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{mpd}}, \code{\link{randomizeMatrix}}
#' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
#' software for the analysis of phylogenetic community structure and trait
#' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
#' @keywords univar
#' @examples
#' data(phylocom)
#' ses.mpd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels")
#' @export ses.mpd
`ses.mpd` <-
  function(
    samp, dis, null.model = c(
      "taxa.labels", "richness", "frequency", "sample.pool",
      "phylogeny.pool", "independentswap", "trialswap"
    ),
    abundance.weighted = FALSE, runs = 999, iterations = 1000) {
    dis <- as.matrix(dis)
    mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted)
    null.model <- match.arg(null.model)
    mpd.rand <- switch(null.model,
                       taxa.labels = t(replicate(runs, mpd(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted))),
                       richness = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
                       frequency = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted))),
                       sample.pool = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
                       phylogeny.pool = t(replicate(runs, mpd(
                         randomizeMatrix(samp, null.model = "richness"),
                         taxaShuffle(dis), abundance.weighted
                       ))),
                       independentswap = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted))),
                       trialswap = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted)))
    )
    mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean) / mpd.rand.sd
    mpd.obs.rank <- apply(
      X = rbind(mpd.obs, mpd.rand), MARGIN = 2,
      FUN = rank
    )[1,]
    mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
    data.frame(
      ntaxa = specnumber(samp), mpd.obs, mpd.rand.mean, mpd.rand.sd, mpd.obs.rank,
      mpd.obs.z, mpd.obs.p = mpd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp)
    )
  }


#' Standardized effect size of MNTD
#' 
#' Standardized effect size of mean nearest taxon distances in communities.
#' When used with a phylogenetic distance matrix, equivalent to -1 times the
#' Nearest Taxon Index (NTI).
#' 
#' The metric used by this function has also been referred to as MNND (mean
#' nearest neighbour distance), and the function was named \code{ses.mnnd} in
#' picante versions < 0.7.
#' 
#' Currently implemented null models (arguments to null.model): \describe{
#' \item{taxa.labels}{ Shuffle distance matrix labels (across all taxa
#' included in distance matrix)} \item{richness}{ Randomize community data
#' matrix abundances within samples (maintains sample species richness)}
#' \item{frequency}{ Randomize community data matrix abundances within species
#' (maintains species occurence frequency)} \item{sample.pool}{ Randomize
#' community data matrix by drawing species from pool of species occurring in
#' at least one community (sample pool) with equal probability}
#' \item{phylogeny.pool}{ Randomize community data matrix by drawing species
#' from pool of species occurring in the distance matrix (phylogeny pool) with
#' equal probability} \item{independentswap}{ Randomize community data matrix
#' with the independent swap algorithm (Gotelli 2000) maintaining species
#' occurrence frequency and sample species richness } \item{trialswap}{
#' Randomize community data matrix with the trial-swap algorithm (Miklos &
#' Podani 2004) maintaining species occurrence frequency and sample species
#' richness } }
#' 
#' @aliases ses.mntd ses.mnnd
#' @param samp Community data matrix
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each
#' species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for
#' independent swap and trial null models)
#' @return A data frame of results for each community \item{ntaxa}{Number of
#' taxa in community} \item{mntd.obs}{Observed MNTD in community}
#' \item{mntd.rand.mean}{Mean MNTD in null communities}
#' \item{mntd.rand.sd}{Standard deviation of MNTD in null communities}
#' \item{mntd.obs.rank}{Rank of observed MNTD vs. null communities}
#' \item{mntd.obs.z}{Standardized effect size of MNTD vs. null communities (=
#' (mntd.obs - mntd.rand.mean) / mntd.rand.sd, equivalent to -NTI)}
#' \item{mntd.obs.p}{P-value (quantile) of observed MNTD vs. null communities
#' (= mntd.obs.rank / runs + 1)} \item{runs}{Number of randomizations}
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{mntd}}, \code{\link{randomizeMatrix}}
#' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
#' software for the analysis of phylogenetic community structure and trait
#' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' ses.mntd(phylocom$sample, cophenetic(phylocom$phylo),null.model="taxa.labels")
#' @export ses.mntd
`ses.mntd` <-
  function(
    samp, dis, null.model = c(
      "taxa.labels", "richness", "frequency", "sample.pool",
      "phylogeny.pool", "independentswap", "trialswap"
    ),
    abundance.weighted = FALSE, runs = 999, iterations = 1000) {
    dis <- as.matrix(dis)
    mntd.obs <- mntd(samp, dis, abundance.weighted)
    null.model <- match.arg(null.model)
    mntd.rand <- switch(null.model,
                        taxa.labels = t(replicate(runs, mntd(samp, taxaShuffle(dis), abundance.weighted))),
                        richness = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
                        frequency = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted))),
                        sample.pool = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
                        phylogeny.pool = t(replicate(runs, mntd(
                          randomizeMatrix(samp, null.model = "richness"),
                          taxaShuffle(dis), abundance.weighted
                        ))),
                        independentswap = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted))),
                        trialswap = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted)))
    )
    mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    mntd.obs.z <- (mntd.obs - mntd.rand.mean) / mntd.rand.sd
    mntd.obs.rank <- apply(
      X = rbind(mntd.obs, mntd.rand), MARGIN = 2,
      FUN = rank
    )[1,]
    mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)
    data.frame(
      ntaxa = specnumber(samp), mntd.obs, mntd.rand.mean, mntd.rand.sd, mntd.obs.rank,
      mntd.obs.z, mntd.obs.p = mntd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp)
    )
  }


psv <- function(samp, tree, compute.var = TRUE, scale.vcv = TRUE) {
  # Make samp matrix a pa matrix
  samp[samp > 0] <- 1

  flag <- 0
  if (is.null(dim(samp))) # if the samp matrix only has one site
  {
    samp <- rbind(samp, samp)
    flag <- 2
  }

  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      tree <- compute.brlen(tree, 1)
    } # If phylo has no given branch lengths
    tree <- prune.sample(samp, tree)
    # Make sure that the species line up
    samp <- samp[, tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix <- vcv.phylo(tree, corr = scale.vcv)
  } else {
    Cmatrix <- tree
    species <- colnames(samp)
    preval <- colSums(samp) / sum(samp)
    species <- species[preval > 0]
    Cmatrix <- Cmatrix[species, species]
    samp <- samp[, colnames(Cmatrix)]
  }

  # numbers of locations and species
  SR <- rowSums(samp)
  nlocations <- dim(samp)[1]
  nspecies <- dim(samp)[2]

  ##################################
  # calculate observed PSVs
  #
  PSVs <- NULL

  for (i in 1:nlocations)
  {
    index <- seq(seq_len(nrow(Cmatrix)))[samp[i,] == 1] # species present
    n <- length(index) # number of species present
    if (n > 1) {
      C <- Cmatrix[index, index] # C for individual locations
      PSV <- (n * sum(diag(as.matrix(C))) - sum(C)) / (n * (n - 1))
    } else {
      PSV <- NA
    }
    PSVs <- c(PSVs, PSV)
  }
  PSVout <- cbind(PSVs, SR)

  if (flag == 2) {
    PSVout <- PSVout[-2,]
    return(PSVout)
  } else {
    if (compute.var == FALSE) {
      return(data.frame(PSVout))
    } else {
      PSVvar <- NULL
      X <- Cmatrix - (sum(sum(Cmatrix - diag(nspecies)))) / (nspecies * (nspecies - 1))
      X <- X - diag(diag(X))

      SS1 <- sum(X * X) / 2

      SS2 <- 0
      for (i in 1:(nspecies - 1)) {
        for (j in (i + 1):nspecies) {
          SS2 <- SS2 + X[i, j] * (sum(X[i,]) - X[i, j])
        }
      }
      SS3 <- -SS1 - SS2

      S1 <- SS1 * 2 / (nspecies * (nspecies - 1))
      S2 <- SS2 * 2 / (nspecies * (nspecies - 1) * (nspecies - 2))

      if (nspecies == 3) {
        S3 <- 0
      } else {
        S3 <- SS3 * 2 / (nspecies *
          (nspecies - 1) *
          (nspecies - 2) *
          (nspecies - 3))
      }

      for (n in 2:nspecies) {
        approxi <- 2 / (n * (n - 1)) * (S1 + (n - 2) * S2 + (n - 2) * (n - 3) * S3)
        PSVvar <- rbind(PSVvar, c(n, approxi))
      }

      vars <- rep(0, nlocations)
      PSVout <- cbind(PSVout, vars)

      for (g in 1:nlocations)
      {
        if (PSVout[g, 2] > 1) {
          PSVout[g, 3] <- PSVvar[PSVout[g, 2] - 1, 2]
        } else {
          PSVout[g, 3] <- NA
        }
      }
      return(data.frame(PSVout))
    }
  }
}

psr <- function(samp, tree, compute.var = TRUE, scale.vcv = TRUE) {
  PSVout <- psv(samp, tree, compute.var, scale.vcv = scale.vcv)
  if (is.null(dim(PSVout)) == TRUE) {
    PSRout <- data.frame(cbind(PSVout[1] * PSVout[2], PSVout[2]))
    names(PSRout) <- c("PSR", "SR")
    rownames(PSRout) <- ""
    return(PSRout)
  } else {
    PSRout <- cbind(PSVout[, 1] * PSVout[, 2], PSVout[, 2])
    if (compute.var != TRUE) {
      colnames(PSRout) <- c("PSR", "SR")
      rownames(PSRout) <- rownames(PSVout)
      return(data.frame(PSRout))
    } else {
      PSRout <- cbind(PSRout, PSVout[, 3] * (PSVout[, 2])^2)
      colnames(PSRout) <- c("PSR", "SR", "vars")
      rownames(PSRout) <- rownames(PSVout)
      return(data.frame(PSRout))
    }
  }
}

pse <- function(samp, tree, scale.vcv = TRUE) {
  flag <- 0
  if (is.null(dim(samp))) # if the samp matrix only has one site
  {
    samp <- rbind(samp, samp)
    flag <- 2
  }

  samp <- as.matrix(samp)

  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      tree <- compute.brlen(tree, 1)
    } # If phylo has no given branch lengths
    tree <- prune.sample(samp, tree)
    # Make sure that the species line up
    samp <- samp[, tree$tip.label, drop = FALSE]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix <- vcv.phylo(tree, corr = scale.vcv)
  } else {
    Cmatrix <- tree
    species <- colnames(samp)
    preval <- colSums(samp) / sum(samp)
    species <- species[preval > 0]
    Cmatrix <- Cmatrix[species, species]
    samp <- samp[, colnames(Cmatrix), drop = FALSE]
  }

  # numbers of locations and species
  SR <- apply(samp > 0, 1, sum)
  nlocations <- dim(samp)[1]

  #################################
  # calculate observed phylogenetic species evenness
  PSEs <- NULL
  for (i in 1:nlocations) {
    index <- seq(1, ncol(Cmatrix))[samp[i,] > 0] # species present
    n <- length(index) # location species richness
    if (n > 1) {
      C <- Cmatrix[index, index] # C for individual locations
      N <- sum(samp[i,]) # location total abundance
      M <- samp[i, samp[i,] > 0] # species abundance column
      mbar <- mean(M) # mean species abundance
      PSE <- (N * t(diag(as.matrix(C))) %*% M - t(M) %*% as.matrix(C) %*% M) / (N^2 - N * mbar) # phylogenetic species evenness
    } else {
      PSE <- NA
    }
    PSEs <- c(PSEs, PSE)
  }
  PSEout <- cbind(PSEs, SR)
  if (flag == 2) {
    PSEout <- PSEout[-2,]
    return(PSEout)
  } else {
    return(data.frame(PSEout))
  }
}


psc <- function(samp, tree, scale.vcv = TRUE) {
  # Make samp matrix a pa matrix
  samp[samp > 0] <- 1
  flag <- 0
  if (is.null(dim(samp))) # if the samp matrix only has one site
  {
    samp <- rbind(samp, samp)
    flag <- 2
  }

  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      tree <- compute.brlen(tree, 1)
    } # If phylo has no given branch lengths
    tree <- prune.sample(samp, tree)
    # Make sure that the species line up
    samp <- samp[, tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix <- vcv.phylo(tree, corr = scale.vcv)
  } else {
    Cmatrix <- tree
    species <- colnames(samp)
    preval <- colSums(samp) / sum(samp)
    species <- species[preval > 0]
    Cmatrix <- Cmatrix[species, species]
    samp <- samp[, colnames(Cmatrix)]
  }

  # numbers of locations and species
  SR <- rowSums(samp)
  nlocations <- dim(samp)[1]

  ##################################
  # calculate observed PSCs
  #
  PSCs <- NULL

  for (i in 1:nlocations)
  {
    index <- seq(seq_len(nrow(Cmatrix)))[samp[i,] == 1] # species present
    n <- length(index) # number of species present
    if (n > 1) {
      C <- Cmatrix[index, index] # C for individual locations
      diag(C) <- -1
      PSC <- 1 - (sum(apply(C, 1, max)) / n)
    } else {
      PSC <- NA
    }
    PSCs <- c(PSCs, PSC)
  }
  PSCout <- cbind(PSCs, SR)
  if (flag == 2) {
    PSCout <- PSCout[-2,]
    return(PSCout)
  } else {
    return(data.frame(PSCout))
  }
}

psv.spp <- function(samp, tree) {
  # Make samp matrix a pa matrix
  samp[samp > 0] <- 1
  if (is.null(dim(samp))) # if the samp matrix only has one site
  {
    samp <- rbind(samp, samp)
  }
  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      tree <- compute.brlen(tree, 1)
    } # If phylo has no given branch lengths
    tree <- prune.sample(samp, tree)
    # Make sure that the species line up
    samp <- samp[, tree$tip.label]
    # Make a correlation matrix of the species pool phylogeny
    Cmatrix <- vcv.phylo(tree, corr = TRUE)
  } else {
    Cmatrix <- tree
    species <- colnames(samp)
    preval <- colSums(samp) / sum(samp)
    species <- species[preval > 0]
    Cmatrix <- Cmatrix[species, species]
    samp <- samp[, colnames(Cmatrix)]
  }
  # reduce given Cmatrix to the species observed in samp
  samp <- samp[rowSums(samp) > 1,] # prune out locations with <2 species

  # cull the species that are not found in the samp set after all communities with 1 species are removed
  preval <- colSums(samp) / sum(samp)
  indexcov <- preval > 0
  Cmatrix <- Cmatrix[indexcov, indexcov]
  samp <- samp[, indexcov]

  obs.PSV <- mean(psv(samp, Cmatrix, compute.var = FALSE)[, 1], na.rm = TRUE)

  # numbers of locations and species
  nspecies <- dim(samp)[2]

  spp.PSVs <- NULL
  for (j in 1:nspecies)
  {
    spp.samp <- samp[, -j]
    spp.Cmatrix <- Cmatrix[-j, -j]
    spp.PSV <- mean(psv(spp.samp, spp.Cmatrix, compute.var = FALSE)[, 1], na.rm = TRUE)
    spp.PSVs <- c(spp.PSVs, spp.PSV)
  }
  spp.PSVout <- (spp.PSVs - obs.PSV) / sum(abs(spp.PSVs - obs.PSV))
  names(spp.PSVout) <- colnames(samp)
  return(spp.PSVout)
}


#' Phylogenetic Species Diversity Metrics
#' 
#' Calculate the bounded phylogenetic biodiversity metrics: phylogenetic
#' species variability, richness, evenness and clustering for one or multiple
#' samples.
#' 
#' \emph{Phylogenetic species variability (PSV)} quantifies how phylogenetic
#' relatedness decreases the variance of a hypothetical unselected/neutral
#' trait shared by all species in a community. The expected value of PSV is
#' statistically independent of species richness, is one when all species in a
#' sample are unrelated (i.e., a star phylogeny) and approaches zero as
#' species become more related. PSV is directly related to mean phylogenetic
#' distance, except except calculated on a scaled phylogenetic covariance
#' matrix. The expected variance around PSV for any sample of a particular
#' species richness can be approximated.  To address how individual species
#' contribute to the mean PSV of a data set, the function \code{psv.spp} gives
#' signed proportions of the total deviation from the mean PSV that occurs
#' when all species are removed from the data set one at a time. The absolute
#' values of these \dQuote{species effects} tend to positively correlate with
#' species prevalence.  \cr \cr \emph{Phylogenetic species richness (PSR)} is
#' the number of species in a sample multiplied by PSV. It can be considered
#' the species richness of a sample after discounting by species relatedness.
#' The value is maximum at the species richness of the sample, and decreases
#' towards zero as relatedness increases. The expected variance around PSR for
#' any sample of a particular species richness can be approximated.  \cr \cr
#' \emph{Phylogenetic species evenness (PSE)} is the metric PSV modified to
#' incorporate relative species abundances. The maximum attainable value of
#' PSE (i.e., 1) occurs only if species abundances are equal and species
#' phylogeny is a star. PSE essentially grafts each individual of a species
#' onto the tip of the phylogeny of its species with branch lengths of zero.
#' \cr \cr \emph{Phylogenetic species clustering (PSC)} is a metric of the
#' branch tip clustering of species across a sample's phylogeny. As PSC
#' increases to 1, species are less related to one another the tips of the
#' phylogeny. PSC is directly related to mean nearest neighbor distance.
#' 
#' @aliases psd psv.spp psv psr pse psc
#' @param samp Community data matrix
#' @param tree A phylo tree object or a phylogenetic covariance matrix
#' @param compute.var Computes the expected variances for PSV and PSR for each
#' community
#' @param scale.vcv Scale the phylogenetic covariance matrix to bound the
#' metric between 0 and 1
#' @return Returns a dataframe of the respective phylogenetic species
#' diversity metric values
#' @note These metrics are bounded either between zero and one (PSV, PSE, PSC)
#' or zero and species richness (PSR); but the metrics asymptotically approach
#' zero as relatedness increases. Zero can be assigned to communities with
#' less than two species, but conclusions drawn from assigning communities
#' zero values need be carefully explored for any data set. The data sets need
#' not be species-community data sets but may be any sample data set with an
#' associated phylogeny.
#' @author Matthew Helmus \email{mrhelmus@@gmail.com}
#' @seealso \code{\link{mpd}} ,\code{\link{mnnd}}, \code{\link{specaccum.psr}}
#' @references Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007)
#' Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' psd(phylocom$sample, phylocom$phylo)
#' 
#' @export psd
psd <- function(samp, tree, compute.var = TRUE, scale.vcv = TRUE) {
  if (is.null(dim(samp))) # if the samp matrix only has one site
  {
    PSDout <- data.frame(c(psv(samp, tree, compute.var, scale.vcv)[1], psc(samp, tree, scale.vcv)[1], psr(samp, tree, compute.var, scale.vcv)[1], pse(samp, tree, scale.vcv)))
    names(PSDout) <- c("PSV", "PSC", "PSR", "PSE", "SR")
    return(PSDout)
  } else {
    if (compute.var == TRUE) {
      PSDout <- cbind(psv(samp, tree, compute.var, scale.vcv)[, c(1, 3)], psc(samp, tree, scale.vcv)[, 1], psr(samp, tree, compute.var, scale.vcv)[, c(1, 3)], pse(samp, tree, scale.vcv))
      colnames(PSDout) <- c("PSV", "var.PSV", "PSC", "PSR", "var.PSR", "PSE", "SR")
    } else {
      PSDout <- cbind(psv(samp, tree, compute.var, scale.vcv)[, 1], psc(samp, tree, scale.vcv)[, 1], psr(samp, tree, compute.var, scale.vcv)[, 1], pse(samp, tree, scale.vcv))
      colnames(PSDout) <- c("PSV", "PSC", "PSR", "PSE", "SR")
    }
    return(data.frame(PSDout))
  }
}


#' Calculate Faith's Phylogenetic Diversity
#' 
#' Calculate the sum of the total phylogenetic branch length for one or
#' multiple samples.
#' 
#' 
#' @param samp Community data matrix
#' @param tree A phylo tree object
#' @param include.root Should the root node be included in all PD calculations
#' (default = TRUE)
#' @return Returns a dataframe of the PD and species richness (SR) values for
#' all samples
#' @note The data sets need not be species-community data sets but may be any
#' sample data set with an associated phylogeny. PD is not statistically
#' independent of species richness, it positively correlates with species
#' richness across samples. The function \code{\link{ses.pd}} compares
#' observed PD to the values expected under various randomizations and allows
#' a way to standardize for unequal richness across samples.
#' 
#' If the root is to be included in all calculations of PD
#' (\code{include.root=TRUE}), the tree must be rooted. Single-species samples
#' will be assigned a PD value equal to the distance from the root to the
#' present.
#' 
#' If the root is not included in all calculations by default
#' (\code{include.root=FALSE}), the tree need not rooted, but in the case of
#' single-species samples the PD will be equal to NA and a warning will be
#' issued.
#' @section Warning : If the root is to be included in all calculations
#' (\code{include.root=TRUE}), the PD of all samples will include the branch
#' length connecting taxa in those samples and the root node of the supplied
#' tree. The root of the supplied tree may not be spanned by any taxa in the
#' sample. If you want the root of your tree to correspond to the most recent
#' ancestor of the taxa actually present in your sample, you should prune the
#' tree before running \code{pd}:
#' 
#' \code{prunedTree <- prune.sample(sample,tree)}
#' @author Matthew Helmus \email{mrhelmus@@gmail.com}, Jonathan Davies
#' \email{davies@@nceas.ucsb.edu}, Steven Kembel
#' \email{steve.kembel@@gmail.com}
#' @seealso \code{\link{psr}}, \code{\link{ses.pd}}
#' @references Faith D.P. (1992) Conservation evaluation and phylogenetic
#' diversity. Biological Conservation, 61, 1-10.
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' pd(phylocom$sample, phylocom$phylo)
#' @export pd
pd <- function(samp, tree, include.root = TRUE) {
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (include.root) {
    # Make sure tree is rooted if needed
    if (!is.rooted(tree)) {
      stop("Rooted tree required to calculate PD with include.root=TRUE argument")
    }
    tree <- node.age(tree)
  }
  species <- colnames(samp)
  SR <- rowSums(ifelse(samp > 0, 1, 0))
  nlocations <- dim(samp)[1]
  PDs <- rep(NA, nlocations)

  for (i in 1:nlocations) {
    present <- species[samp[i,] > 0]
    treeabsent <- tree$tip.label[which(!(tree$tip.label %in% present))]
    # Check that sample has species; If no species are present, PD = 0
    if (length(present) == 0) {
      PDs[i] <- 0
    }
      # Check if there is only ONE species in sample
    else if (length(present) == 1) {
      # If tree is not rooted, parse error message
      if (!is.rooted(tree) || !include.root) {
        warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.")
        PDs[i] <- NA
      }
        # Else the PD is the node age of that single species present
      else {
        PDs[i] <- tree$ages[which(tree$edge[, 2] ==
                                    which(tree$tip.label == present))]
      }
    }
      # If there are no absent species (all tips are in sample) then PD
      # is the sum of all edges in the tree
    else if (length(treeabsent) == 0) {
      PDs[i] <- sum(tree$edge.length)
    }
      # Otherwise, we will need to remove absent species
    else {
      # Make a SubTree with only present species
      sub.tree <- drop.tip(tree, treeabsent)
      # If the tree is rooted, you need to account for different in
      # maximum branch length between subtree and original
      if (include.root) {
        # Make sure tree is rooted if needed
        if (!is.rooted(tree)) {
          stop("Rooted tree required to calculate PD with include.root=TRUE argument")
        }
        # Calculate diff btw max depth original and max depth subtree
        sub.tree.depth <- max(node.age(sub.tree)$ages)
        orig.tree.depth <- max(tree$ages[which(tree$edge[, 2] %in% which(tree$tip.label %in% present))])
        # PD is the sum of edges in subtree, plus any difference in the
        # max distance (of only present tips) between original and sub
        PDs[i] <- sum(sub.tree$edge.length) + (orig.tree.depth -
          sub.tree.depth)
      }
        # If root is not included, just add all edges
      else {
        PDs[i] <- sum(sub.tree$edge.length)
      }
    }
  }
  PDout <- data.frame(PD = PDs, SR = SR)
  rownames(PDout) <- rownames(samp)
  return(PDout)
}


#' Standardized effect size of PD
#' 
#' Standardized effect size of phylogenetic diversity (Faith's PD) in
#' communities.
#' 
#' Currently implemented null models (arguments to null.model): \describe{
#' \item{taxa.labels}{ Shuffle taxa labels across tips of phylogeny (across
#' all taxa included in phylogeny)} \item{richness}{ Randomize community data
#' matrix abundances within samples (maintains sample species richness)}
#' \item{frequency}{ Randomize community data matrix abundances within species
#' (maintains species occurence frequency)} \item{sample.pool}{ Randomize
#' community data matrix by drawing species from pool of species occurring in
#' at least one community (sample pool) with equal probability}
#' \item{phylogeny.pool}{ Randomize community data matrix by drawing species
#' from pool of species occurring in the phylogeny (phylogeny pool) with equal
#' probability} \item{independentswap}{ Randomize community data matrix with
#' the independent swap algorithm (Gotelli 2000) maintaining species
#' occurrence frequency and sample species richness } \item{trialswap}{
#' Randomize community data matrix with the trial-swap algorithm (Miklos &
#' Podani 2004) maintaining species occurrence frequency and sample species
#' richness } }
#' 
#' @param samp Community data matrix
#' @param tree Phylogeny (phylo object)
#' @param null.model Null model to use (see Details section for description)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for
#' independent swap and trial null models)
#' @param include.root Include distance to root node in calculation of PD (see
#' documentation in \code{\link{pd}} function
#' @return A data frame of results for each community \item{ntaxa}{Number of
#' taxa in community} \item{pd.obs}{Observed PD in community}
#' \item{pd.rand.mean}{Mean PD in null communities} \item{pd.rand.sd}{Standard
#' deviation of PD in null communities} \item{pd.obs.rank}{Rank of observed PD
#' vs. null communities} \item{pd.obs.z}{Standardized effect size of PD vs.
#' null communities (= (pd.obs - pd.rand.mean) / pd.rand.sd)}
#' \item{pd.obs.p}{P-value (quantile) of observed PD vs. null communities (=
#' mpd.obs.rank / runs + 1)} \item{runs}{Number of randomizations}
#' @author Steven Kembel <steve.kembel@@gmail.com>
#' @seealso \code{\link{pd}}, \code{\link{randomizeMatrix}}
#' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
#' software for the analysis of phylogenetic community structure and trait
#' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
#' 
#' Proches, S., Wilson, J.R.U. and Cowling, R.M. 2006. How much evolutionary
#' history in a 10 x 10m plot? Proceedings of Royal Society of London B,
#' Biological Sciences 273:1143-1148.
#' @keywords univar
#' @examples
#' 
#' data(phylocom)
#' ses.pd(phylocom$sample, phylocom$phylo, null.model="taxa.labels", runs=99)
#' @export ses.pd
ses.pd <- function(
  samp, tree, null.model = c(
    "taxa.labels", "richness", "frequency",
    "sample.pool", "phylogeny.pool", "independentswap", "trialswap"
  ),
  runs = 999, iterations = 1000, include.root = TRUE) {
  if (include.root == TRUE) {
    pd.obs <- as.vector(pd(samp, tree, include.root = TRUE)$PD)
    null.model <- match.arg(null.model)
    pd.rand <- switch(null.model,
                      taxa.labels = t(replicate(runs, as.vector(pd(samp, tipShuffle(tree), include.root = TRUE)$PD))),
                      richness = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "richness"), tree, include.root = TRUE)$PD))),
                      frequency = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "frequency"), tree, include.root = TRUE)$PD))),
                      sample.pool = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "richness"), tree, include.root = TRUE)$PD))),
                      phylogeny.pool = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "richness"),
                                                                      tipShuffle(tree),
                                                                      include.root = TRUE
                      )$PD))),
                      independentswap = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "independentswap", iterations), tree, include.root = TRUE)$PD))),
                      trialswap = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "trialswap", iterations), tree, include.root = TRUE)$PD)))
    )
    pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    pd.obs.z <- (pd.obs - pd.rand.mean) / pd.rand.sd
    pd.obs.rank <- apply(
      X = rbind(pd.obs, pd.rand), MARGIN = 2,
      FUN = rank
    )[1,]
    pd.obs.rank <- ifelse(is.na(pd.rand.mean), NA, pd.obs.rank)
    return(data.frame(
      ntaxa = specnumber(samp), pd.obs, pd.rand.mean, pd.rand.sd, pd.obs.rank,
      pd.obs.z, pd.obs.p = pd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp)
    ))
  }

  if (include.root == FALSE) {
    pd.obs <- as.vector(pd(samp, tree, include.root = FALSE)$PD)
    null.model <- match.arg(null.model)
    pd.rand <- switch(null.model,
                      taxa.labels = t(replicate(runs, as.vector(pd(samp, tipShuffle(tree), include.root = FALSE)$PD))),
                      richness = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "richness"), tree, include.root = FALSE)$PD))),
                      frequency = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "frequency"), tree, include.root = FALSE)$PD))),
                      sample.pool = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "richness"), tree, include.root = FALSE)$PD))),
                      phylogeny.pool = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "richness"),
                                                                      tipShuffle(tree),
                                                                      include.root = FALSE
                      )$PD))),
                      independentswap = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "independentswap", iterations), tree, include.root = FALSE)$PD))),
                      trialswap = t(replicate(runs, as.vector(pd(randomizeMatrix(samp, null.model = "trialswap", iterations), tree, include.root = FALSE)$PD)))
    )
    pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    pd.obs.z <- (pd.obs - pd.rand.mean) / pd.rand.sd
    pd.obs.rank <- apply(
      X = rbind(pd.obs, pd.rand), MARGIN = 2,
      FUN = rank
    )[1,]
    pd.obs.rank <- ifelse(is.na(pd.rand.mean), NA, pd.obs.rank)
    return(data.frame(
      ntaxa = specnumber(samp), pd.obs, pd.rand.mean, pd.rand.sd, pd.obs.rank,
      pd.obs.z, pd.obs.p = pd.obs.rank / (runs + 1), runs = runs, row.names = row.names(samp)
    ))
  }
}
