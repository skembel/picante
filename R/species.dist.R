##' Species co-occurrence distances
##' 
##' Compute interspecific distances based on patterns of species co-occurrence
##' in communities.
##' 
##' Currently implemented co-occurrence measures (arguments to metric):
##' \describe{ \item{cij}{ Schoener's index of co-occurrence } \item{jaccard}{
##' Jaccard index of co-occurrence } \item{checkerboard}{ Checkerboard index of
##' co-occurrence } \item{doij}{ DOij index of co-occurrence } }
##' 
##' @param x Community data matrix
##' @param metric Co-occurrence metric to use (see Details section for
##' description)
##' @return A \code{dist} object with co-occurrences among all species pairs
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @seealso \code{\link[vegan]{vegdist}}
##' @references Hardy, O.J. 2008. Testing the spatial phylogenetic structure of
##' local communities: statistical performances of different null models and
##' test statistics on a locally neutral community. Journal of Ecology
##' 96:914-926.
##' @keywords univar
##' @export species.dist
`species.dist` <-
  function(x, metric = c("cij", "jaccard", "checkerboard", "doij")) {
    metric <- match.arg(metric)
    if (identical(metric, "checkerboard")) {
      # Gotelli 2000: Checker = Sum (Si - Q)(Sk - Q) / ((R*(R-1))/2)
      # where Si = total for row(species) i, R = num rows(spp), Q = num sites where both spp present
      x <- decostand(x, method = "pa")
      Nsites <- dim(x)[1]
      S <- apply(x, 2, sum)
      R <- length(S)
      Checker.ij <- matrix(nrow = R, ncol = R, dimnames = list(colnames(x), colnames(x)))
      for (i in 1:R) {
        for (j in 1:R) {
          Q <- sum(x[, i] * x[, j])
          Checker.ij[i, j] <- ((S[i] - Q) * (S[j] - Q)) / ((R * (R - 1)) / 2)
        }
      }
      return(as.dist(Checker.ij))
    }
    if (identical(metric, "cij")) {
      # Schoener index of co-occurrence
      x <- decostand(x, method = "total", MARGIN = 2)
      return(1 - (0.5 * dist(t(x), method = "manhattan")))
    }
    if (identical(metric, "jaccard")) {
      return(1 - vegdist(t(sortColumns(x)), method = "jaccard"))
    }
    if (identical(metric, "doij")) {
      # Hardy's standardized version of checkerboard
      # doij = (Pij - Pi*Pj)/(Pi*Pj)
      x <- as.matrix(decostand(x, method = "pa"))
      Nsites <- dim(x)[1]
      P <- apply(x, 2, sum) / Nsites
      N <- length(P)
      doij <- matrix(nrow = N, ncol = N, dimnames = list(colnames(x), colnames(x)))
      for (i in 1:N - 1) {
        for (j in (i + 1):N) {
          Pij <- sum(x[, i] * x[, j]) / Nsites
          doij[i, j] <- ((Pij - (P[i] * P[j])) / (P[i] * P[j]))
        }
      }
      return(as.dist(t(doij)))
    }
  }
