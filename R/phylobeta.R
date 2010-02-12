#phylobeta diversity functions

#comdist
#mean pairwise distance among taxa from two communities
comdist <- function(comm, dis, abundance.weighted=FALSE) {

    x <- as.matrix(comm)
    dat <- match.comm.dist(comm, dis)
    x <- dat$comm
    dis <- as.matrix(dat$dist)
    if (!abundance.weighted) {
        x <- decostand(x, method="pa")
    }
    N <- dim(x)[1]
    S <- dim(x)[2]
    x <- decostand(x, method="total", MARGIN=1)

	comdist <- matrix(nrow=N,ncol=N)    
	for (l in 1:(N-1)) {
		for (k in 2:N) {
			comdist[k,l] <-
			sum ( dis * outer(as.vector(t(x[k,])),as.vector(t(x[l,]))) )
		}
	}

    row.names(comdist) <- row.names(x)
    colnames(comdist) <- row.names(x)
    return(as.dist(comdist))

}


#comdistnt
#mean distance to closest relative between taxa from two communities
comdistnt <- function(comm, dis, abundance.weighted=FALSE, exclude.conspecifics=FALSE) {

    dat <- match.comm.dist(comm, dis)
    comm <- dat$comm
    dis <- dat$dist
    N <- dim(comm)[1]
    comm <- decostand(comm, method="total", MARGIN=1)    
    comdisnt <- matrix(nrow=N,ncol=N)
    for (i in 1:(N-1)) {
        for (j in (i+1):N) {
            sppInSample1 <- colnames(comm[i, comm[i, ] > 0, drop=FALSE])
            sppInSample2 <- colnames(comm[j, comm[j, ] > 0, drop=FALSE])
            if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 1)) {
                sample.dis <- dis[sppInSample1, sppInSample2, drop=FALSE]
                if (exclude.conspecifics) {
                    sample.dis[sample.dis==0] <- NA
                }
                #TODO fix min throws errors on empty set
                sample1NT <- apply(sample.dis,1,min,na.rm=TRUE)
                sample1NT[sample1NT == Inf] <- NA
                sample2NT <- apply(sample.dis,2,min,na.rm=TRUE)
                sample2NT[sample2NT == Inf] <- NA     
                if (abundance.weighted) {
                    sample1.weights <- as.numeric(comm[i,sppInSample1])
                    sample2.weights <- as.numeric(comm[j,sppInSample2])
                    if (any(is.na(sample1NT))) {
                        miss <- which(is.na(sample1NT))
                        sample1NT <- sample1NT[-miss]
                        sample1.weights <- sample1.weights[-miss]
                        sample1.weights <- sample1.weights / sum(sample1.weights)
                    }
                    if (any(is.na(sample2NT))) {
                        miss <- which(is.na(sample2NT))
                        sample2NT <- sample2NT[-miss]
                        sample2.weights <- sample2.weights[-miss]
                        sample2.weights <- sample2.weights / sum(sample2.weights)
                    }
                    sampleNT <- c(sample1NT, sample2NT)
                    sample.weights <- c(sample1.weights,sample2.weights)
                    comdisnt[i,j] <- weighted.mean(sampleNT, sample.weights, na.rm=TRUE)
                }
                else
                {
                    comdisnt[i,j] <- mean(c(sample1NT,sample2NT), na.rm=TRUE)                
                }
            }
            else{
                comdisnt[i,j] <- NA
            }
        }
    }
    rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
    return(as.dist(t(comdisnt)))
}

#unifrac
unifrac <- function (comm, tree) 
{
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }
    
    comm <- as.matrix(comm)
    s <- nrow(comm)
    phylodist <- matrix(NA, s, s)
    rownames(phylodist) <- rownames(comm)
    colnames(phylodist) <- rownames(comm)
        
    comm_comb<-matrix(NA,s*(s-1)/2,ncol(comm))
    colnames(comm_comb)<-colnames(comm)
    
    i<-1
    for (l in 1:(s - 1))
    {
    	for (k in (l + 1):s)
    	{comm_comb[i,]<-comm[l, ] + comm[k, ]
    	i<-i+1}
    }
        	
    pdcomm<-pd(comm, tree)
	pdcomm_comb<-pd(comm_comb, tree)
    
    i<-1
    for (l in 1:(s - 1)) {
        pdl <- pdcomm[l,"PD"]
        for (k in (l + 1):s) {
            pdk <- pdcomm[k,"PD"]
            pdcomb <- pdcomm_comb[i,"PD"]
            pdsharedlk <- pdl + pdk - pdcomb
            phylodist[k, l] = (pdcomb - pdsharedlk) / pdcomb
            i<-i+1
        }
    }
    return(as.dist(phylodist))
}
