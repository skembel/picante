#' @title match.phylo...
#' @aliases match.phylo.data match.phylo.comm match.phylo.dist match.comm.dist
#' @description
#' These functions compare taxa present in phylogenies with community or
#' trait data sets, pruning and sorting the two kinds of data to match one 
#' another for subsequent analysis.
#' 
#' @usage 
#' 
#' match.phylo.comm(phy, comm)
#' match.phylo.data(phy, data)
#' match.comm.dist(comm, dist)
#' 
#' @param phy A phylogeny object of class phylo
#' @param comm Community data matrix
#' @param data A data object - a vector (with names matching phy) or a
#' data.frame or matrix (with row names matching phy)
#' @param dist A distance matrix - a dist or matrix object
#' 
#' @details
#' A common pitfall in comparative analyses in R is that taxa labels
#' are assumed to match between phylogenetic and other data sets.
#' These functions prune a phylogeny and community or trait data set
#' to match one another, reporting taxa that are missing from one data
#' set or the other.Taxa names for phylogeny objects are taken from
#' the phylogeny's tip labels. Taxa names for community data are taken
#' from the column names. Taxa names for trait data are taken from the element
#' names (vector) or row names (data.frame or matrix).
#' Taxa names for distance data are taken from column/row names of
#' the distance matrix/dist object. If community data lack taxa names, 
#' the function will issue a warning
#' and no result will be returned, since the community-phylogenetic
#' analyses in \code{picante} require named taxa in the community data set. 
#' If trait data or distance matrix lack names, a warning is issued and
#' the data are assumed to be sorted in the same order as the phylogeny's tip 
#' labels or community's column labels.These utility functions are used
#' by several functions that assume taxa labels
#' in phylogeny and data match, including \code{\link{Kcalc}},
#' \code{\link{phylosignal}}, and \code{\link{raoD}}.
#'
#' @author{ Steven Kembel <steve.kembel@gmail.com> }
#' @seealso \code{\link{prune.missing}}, \code{\link{prune.sample}}
#' 
#' @examples
#' data(phylocom)
#' match.phylo.comm(phylocom$phylo, phylocom$sample)
#' match.phylo.data(phylocom$phylo, phylocom$traits[1:10,])
#' @keywords univar
#' @export match.phylo.comm
#' @export match.phylo.data


match.phylo.comm <- function(phy, comm) {
  if (!(is.data.frame(comm) || is.matrix(comm))) {
    stop("Community data should be a data.frame or matrix with
     samples in rows and taxa in columns")
  }


  res <- list()
  phytaxa <- phy$tip.label
  commtaxa <- colnames(comm)

  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names,
    these are required to match phylogeny and community data")
  }

  if (!all(commtaxa %in% phytaxa)) {
    print("Dropping taxa from the community because
     they are not present in the phylogeny:")
    print(setdiff(commtaxa, phytaxa))
    comm <- comm[, intersect(commtaxa, phytaxa)]
    commtaxa <- colnames(comm)
  }

  if (any(!(phytaxa %in% commtaxa))) {
    print("Dropping tips from the tree
     because they are not present in the community data:")
    print(setdiff(phytaxa, commtaxa))
    res$phy <- prune.sample(comm, phy)
  } else {
    res$phy <- phy
  }

  res$comm <- comm[, res$phy$tip.label]
  return(res)
}

match.phylo.data <- function(phy, data) {
  res <- list()
  phytaxa <- phy$tip.label

  if (is.vector(data)) {
    datataxa <- names(data)
    if (is.null(datataxa)) {
      warning("Data set lacks taxa names, these are
       required to match phylogeny and data. Data are returned unsorted.
        Assuming that data and phy$tip.label are in the same order!")
      return(list(phy = phy, data = data))
    }

    if (!all(datataxa %in% phytaxa)) {
      print("Dropping taxa from the data because
       they are not present in the phylogeny:")
      print(setdiff(datataxa, phytaxa))
      data <- data[intersect(datataxa, phytaxa)]
      datataxa <- names(data)
    }

    if (any(!(phytaxa %in% datataxa))) {
      print("Dropping tips from the tree because
       they are not present in the data:")
      print(setdiff(phytaxa, datataxa))
      res$phy <- drop.tip(phy, setdiff(phytaxa, datataxa))
    } else {
      res$phy <- phy
    }

    res$data <- data[res$phy$tip.label]
    return(res)

    ## if they are not vector

  } else if (is.data.frame(data) || is.matrix(data)) {
    dataclass <- class(data)
    data <- as.matrix(data)
    datataxa <- rownames(data)

    if (is.null(datataxa)) {
      warning("Data set lacks taxa (row) names,
       these are required to match phylogeny and data.
        Data are returned unsorted. Assuming that data rows
         and phy$tip.label are in the same order!")
      return(list(phy = phy, data = data))
    }

    if (!all(datataxa %in% phytaxa)) {
      print("Dropping taxa from the data because
       they are not present in the phylogeny:")
      print(setdiff(datataxa, phytaxa))
      data <- data[intersect(datataxa, phytaxa)]
      datataxa <- rownames(data)

    }

    if (any(!(phytaxa %in% datataxa))) {
      print("Dropping tips from the tree because
       they are not present in the data:")
      print(setdiff(phytaxa, datataxa))
      res$phy <- drop.tip(phy, setdiff(phytaxa, datataxa))
    } else {
      res$phy <- phy
    }

    if (dataclass == "data.frame") {

      res$data <- data[res$phy$tip.label,]

    } else {
      res$data <- data[res$phy$tip.label,]
    }
    return(res)
  } else {
    stop("Data must be a vector, data.frame, or matrix")
  }
}


match.comm.dist <- function(comm, dist) {
  res <- list()

  commtaxa <- colnames(comm)

  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names,
     these are required to match distance matrix and community data")
  }

  disclass <- dist
  dist <- as.matrix(dist)

  distaxa <- rownames(dist)

  if (is.null(distaxa)) {
    warning("Distance matrix lacks taxa names, these are required
     to match community and distance matrix. Data are returned unsorted.
      Assuming that distance matrix and community data taxa columns
       are in the same order!")
    if (inherits(disclass, "dist")) {
      return(list(comm = comm, dist = as.dist(dist)))
    } else {
      return(list(comm = comm, dist = dist))
    }
  }

  if (!all(distaxa %in% commtaxa)) {
    print("Dropping taxa from the distance matrix because
     they are not present in the community data:")
    print(setdiff(distaxa, commtaxa))
    dist <- dist[intersect(distaxa, commtaxa), intersect(distaxa, commtaxa)]
    distaxa <- rownames(dist)
  }

  if (any(!(commtaxa %in% distaxa))) {
    print("Dropping taxa from the community because
     they are not present in the distance matrix:")
    print(setdiff(commtaxa, distaxa))
    res$comm <- comm[, intersect(commtaxa, distaxa)]
  } else {
    res$comm <- comm
  }

  if (inherits(disclass, "dist")) {
    res$dist <- as.dist(dist[colnames(comm), colnames(comm)])
  } else {
    res$dist <- dist[colnames(comm), colnames(comm)]
  }
  return(res)
}
