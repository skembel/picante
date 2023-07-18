##' Null models for community data matrix randomization
##' 
##' Various null models for randomizing community data matrices
##' 
##' Currently implemented null models (arguments to null.model): \describe{
##' \item{frequency}{ Randomize community data matrix abundances within species
##' (maintains species occurence frequency)} \item{richness}{ Randomize
##' community data matrix abundances within samples (maintains sample species
##' richness)} \item{independentswap}{ Randomize community data matrix with the
##' independent swap algorithm (Gotelli 2000) maintaining species occurrence
##' frequency and sample species richness } \item{trialswap}{ Randomize
##' community data matrix with the trial-swap algorithm (Miklos & Podani 2004)
##' maintaining species occurrence frequency and sample species richness } }
##' 
##' @param samp Community data matrix
##' @param null.model Null model to use (see Details section for description)
##' @param iterations Number of independent or trial-swaps to perform
##' @return Randomized community data matrix
##' @author Steven Kembel <steve.kembel@@gmail.com>
##' @references Gotelli, N.J. 2000. Null model analysis of species
##' co-occurrence patterns. Ecology 81: 2606-2621
##' 
##' Miklos I. & Podani J. 2004. Randomization of presence-absence matrices:
##' Comments and new algorithms. Ecology 85: 86-92.
##' @keywords manip
##' @examples
##' data(phylocom)
##' randomizeMatrix(phylocom$sample, null.model="richness")
##' @export randomizeMatrix
##' @useDynLib picante, .registration = TRUE

`randomizeMatrix` <-
  function(samp, null.model=c("frequency","richness","independentswap","trialswap"),
           iterations=1000)
  {
    
    samp <- as.matrix(samp)
    null.model <- match.arg(null.model)
    #for independent and trial swap - how to swap abundances?
    #abundance.swap=c("species","sites","both")
    #abundance swap = 0 (within species), 1 (within sites), 2 (random)    
    #abundance.swap <- match.arg(abundance.swap)
    #abundance.swap <- match(abundance.swap,c("species","sites","random"))
    
    if (identical(null.model,"frequency")) {
      ret1 <- .C("frequency", m=as.numeric(samp), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
      return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))
    }
    
    if (identical(null.model,"richness")) {
      ret1 <- .C("richness", m=as.numeric(samp), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
      return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))
    }
    
    if (identical(null.model,"independentswap")) 
    {
      ret1 <- .C("independentswap", m=as.numeric(samp), as.integer(iterations), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
      return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))        
    }
    
    if (identical(null.model,"trialswap")) 
    {
      ret1 <- .C("trialswap", m=as.numeric(samp), as.integer(iterations), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
      return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))        
    }
    
  }

