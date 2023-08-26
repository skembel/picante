##' Convert Phylocom sample to community data matrix
##' 
##' Convert a Phylocom database-format sample to community data matrix.
##' 
##' @param x Phylocom sample formatted data frame, a data frame with three
##' 
##' |  Community name  |  Species abundance  | Species name |
##' | --- | --- | --- |
##' | 0 | 0 |   0  |
##' 
##'
##' @author Steven Kembel <steve.kembel@@gmail.com> and Cam Webb
##' <cwebb@@oeb.harvard.edu>
##' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
##' @keywords IO
##' @export sample2matrix
`sample2matrix` <-
  function(x) {
    # turns a phylocom-format sample file into a sample X species matrix
    # assumes a phylocom-format sample:
    # no header line data
    # column 1 = plot, column 2 = abundance, column 3 = species
    # to load a phylocom-format file directly, try:
    # sample2matrix(read.delim2(file="FILENAME",header=F))
    colnames(x) <- c("plot", "abund", "id")
    y <- tapply(x$abund, list(x$plot, x$id), sum)
    y[is.na(y)] <- 0
    as.data.frame(y)
  }



##' Convert community data matrix to Phylocom sample
##' 
##' Converts a community data matrix to a Phylocom database-format community
##' sample
##' 
##' 
##' @param z Community data matrix
##' @return Phylocom database-format community sample
##' @author Steven Kembel <steve.kembel@@gmail.com> and Cam Webb
##' <cwebb@@oeb.harvard.edu>
##' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
##' @keywords manip
##' @examples
##' 
##' data(phylocom)
##' matrix2sample(phylocom$sample)
##' 
##' @export matrix2sample
`matrix2sample` <-
  function(z) {
    temp <- data.frame(
      expand.grid(dimnames(provideDimnames(z)))[1:2],
      as.vector(as.matrix(z))
    )
    temp <- temp[(temp[, 3] > 0) & !is.na(temp[, 3]), ]
    temp <- temp[sort.list(temp[, 1]), ]
    data.frame(plot = temp[, 1], abund = temp[, 3], id = temp[, 2])
  }




##' Read Phylocom sample
##' 
##' Reads a Phylocom sample file and converts to a community data matrix
##' 
##' 
##' @param filename Phylocom sample file path
##' @return Community data matrix
##' @author Steven Kembel <skembel> and Cam Webb <cwebb@@oeb.harvard.edu>
##' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
##' @keywords IO
##' @export readsample
`readsample` <-
  function(filename = "") {
    x <- read.table(file = filename, header = FALSE, sep = "\t", col.names = c(
      "plot",
      "abund", "id"
    ))
    sample2matrix(x)
  }




##' Write a Phylocom community sample file
##' 
##' Write a community data matrix to a Phylocom community sample file
##' 
##' 
##' @param community Community data matrix
##' @param filename Filename path
##' @author Steven Kembel <steve.kembel@@gmail.com> and Cam Webb
##' <cwebb@@oeb.harvard.edu>
##' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
##' @keywords file
##' @export writesample
`writesample` <-
  function(community, filename = "") {
    write.table(matrix2sample(community),
      file = filename, append = FALSE,
      sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }




##' Write a Phylocom traits formatted file
##' 
##' Write a Phylocom traits formatted file
##' 
##' 
##' @param trt Data frame containing trait data
##' @param file Filename path
##' @param bin Vector index of trait columns to be treated as binary
##' @param sigd Significant digits for output
##' @author David Ackerly <dackerly@@berkeley.edu> and Steven Kembel
##' <steve.kembel@@gmail.com>
##' @references Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:
##' software for the analysis of phylogenetic community structure and trait
##' evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}.
##' @keywords file
##' @export writetraits
`writetraits` <-
  function(trt, file = "", bin = NULL, sigd = 3) {
    head <- matrix("3", ncol = length(names(trt)), nrow = 1)
    if (!is.null(bin)) {
      head[bin] <- 0
    }
    write.table(data.frame("type", head),
      file = file, sep = "\t",
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(matrix(c("name", colnames(trt)), nrow = 1),
      file = file, sep = "\t", append = TRUE, quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )
    write.table(signif(trt, sigd),
      sep = "\t", file = file, append = TRUE, quote = FALSE,
      row.names = TRUE, col.names = FALSE
    )
  }

#' @title  phylocom
#' @name phylocom
#' @description one row with NA
#' @docType data
#' @usage data(phylocom)
#' @keywords datasets
NULL
