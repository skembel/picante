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
