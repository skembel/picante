

#' @title Net relatedness index (NRI)
#' @description
#' Metric based on the mean phylogenetic distance (MPD) to measures the standardized effect size of the MPD
#' @describeIn ses.mpd 
#' @export NRI
#' @examples
#' data(phylocom)
#' NRI(samp = phylocom$sample, dis = cophenetic(phylocom$phylo),null.model="taxa.labels")
#' NTI(samp = phylocom$sample, dis = cophenetic(phylocom$phylo),null.model="taxa.labels")

#'
#' @param samp Community data matrix
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each
#' species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for
#' independent swap and trial null models)
NRI <- function(samp,
                dis,
                null.model = c(
                  "taxa.labels",
                  "richness",
                  "frequency",
                  "sample.pool",
                  "phylogeny.pool",
                  "independentswap",
                  "trialswap"
                ),
                abundance.weighted = FALSE,
                runs = 999,
                iterations = 1000) {
  mpd.result <-
    picante::ses.mpd(samp, dis, null.model, abundance.weighted, runs, iterations)
  
  NRI.result <- as.matrix(-1 *
                            ((mpd.result[, 2] - mpd.result[, 3]) /
                               mpd.result[, 4]))
  rownames(NRI.result) <- row.names(mpd.result)
  colnames(NRI.result) <- 'NRI'
  as.data.frame(NRI.result)
  NRI.result
  
}



#' @title nearest taxon index (NTI) 
#' @param samp Community data matrix
#' @param dis Distance matrix (generally a phylogenetic distance matrix)
#' @param null.model Null model to use (see Details section for description)
#' @param abundance.weighted Should mean nearest taxon distances for each
#' species be weighted by species abundance? (default = FALSE)
#' @param runs Number of randomizations
#' @param iterations Number of iterations to use for each randomization (for
#' @export NTI
#' @description
#' calculates the mean nearest phylogenetic neighbor among
#'  the individuals in a community
#' @references
#' {Chai, Y., Yue, M., Liu, X. et al. Patterns of taxonomic, 
#' phylogenetic diversity during a long-term succession of forest
#'  on the Loess Plateau, China: insights into assembly process. 
#'  Sci Rep 6, 27087 (2016). https://doi.org/10.1038/srep27087}
#'
#' @describeIn NRI 

NTI <- function(samp,
                dis,
                null.model = c(
                  "taxa.labels",
                  "richness",
                  "frequency",
                  "sample.pool",
                  "phylogeny.pool",
                  "independentswap",
                  "trialswap"
                ),
                abundance.weighted = FALSE,
                runs = 999,
                iterations = 1000) {
  mntd.result <-
    picante::ses.mntd(samp, dis, null.model, abundance.weighted, runs, iterations)
  
  NTI.result <- as.matrix(-1 *
                            ((mntd.result[, 2] - mntd.result[, 3]) /
                               mntd.result[, 4]))
  rownames(NTI.result) <- row.names(mntd.result)
  colnames(NTI.result) <- 'NTI'
  NTI.result
}
