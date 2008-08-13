.First.lib <- function(lib,pkg) {
  library.dynam("picante",pkg,lib)
} 

`randomizeSample` <-
function(samp, null.model=c("frequency","richness","independentswap","trialswap"), lang=c("C","R"), it=1000) {

    samp <- as.matrix(samp)
    lang <- match.arg(lang)
    null.model <- match.arg(null.model)
    
	if (identical(null.model,"frequency")) {
        if (identical(lang,"R")) {
    	    return(data.frame(apply(samp,2,sample),row.names=row.names(samp)))
    	} else {
    	    ret1 <- .C("frequency", m=as.numeric(samp), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
    	    return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))
    	}
	}

	if (identical(null.model,"richness")) {
	    if (identical(lang,"R")) {
            return(t(data.frame(apply(samp,1,sample),row.names=colnames(samp))))
       	} else {
    	    ret1 <- .C("richness", m=as.numeric(samp), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
    	    return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))
    	}
	}

	if (identical(null.model,"independentswap")) 
    {
        if (identical(lang,"R")) {
            #check for presence-absence and warn until abundance implemented
            x <- decostand(samp, "pa")
            if (!identical(x,samp)) stop("Null model currently requires a presence-absence matrix.")
            sppFreq <- apply(x, 2, sum)/ncol(x)
            siteRichness <- apply(x, 1, sum)
            sampleList <- vector()
            sppList <- vector()
            for (siteNum in 1:length(siteRichness)) {
                sampleList <- c(sampleList, rep(names(siteRichness)[siteNum], 
                    siteRichness[siteNum]))
                sppList <- c(sppList, sample(names(sppFreq), siteRichness[siteNum], 
                    replace = FALSE, prob = sppFreq))
            }
            shuffledList <- data.frame(sample = sampleList, species = sppList, 
                p = rep(1, length(sppList)))
            shuffledMatrix <- tapply(shuffledList$p, list(shuffledList$sample, 
                shuffledList$species), sum)
            shuffledMatrix[is.na(shuffledMatrix)] <- 0
            if (identical(dim(shuffledMatrix),dim(x))) {
                return(shuffledMatrix[rownames(x), ])        
            } else {
                mergedFrame <- merge(shuffledMatrix, x, all.y = TRUE)
                mergedFrame[is.na(mergedFrame)] <- 0
                return(mergedFrame[rownames(x),colnames(x)])
            }        
        } else {
            ret1 <- .C("independentswap", m=as.numeric(samp), as.integer(it), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
            return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))        
        }
	}
	
    if (identical(null.model,"trialswap")) 
    {
        if (identical(lang,"R")) {
            stop("Trial swap implemented in C only")
        } else {
            ret1 <- .C("trialswap", m=as.numeric(samp), as.integer(it), as.integer(nrow(samp)), as.integer(ncol(samp)), PACKAGE="picante")
            return(matrix(ret1$m,nrow=nrow(samp),dimnames=list(rownames(samp),colnames(samp))))        
        }
	}

}
