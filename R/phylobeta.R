#phylobeta diversity functions

#comdist
#mean pairwise distance among taxa from two communities
comdist <- function(comm,dis,abundance.weighted=FALSE,exclude.conspecifics=FALSE) {
    N <- dim(comm)[1]
    #relativize abundances for inter-sample comparisons
    comm <- decostand(comm,method="total",MARGIN=1)
    comdis <- matrix(nrow=N,ncol=N)
    for (i in 1:(N-1)) {
        for (j in (i+1):N) {
            sppInSample1 <- names(comm[i, comm[i, ] > 0])
            sppInSample2 <- names(comm[j, comm[j, ] > 0])
            if ((length(sppInSample1) > 1) && (length(sppInSample2) > 1)) {
                sample.dis <- dis[sppInSample1, sppInSample2]
                if (exclude.conspecifics) {
                    sample.dis[sample.dis == 0] <- NA
                }                
                if (abundance.weighted) {
                    sample.weights <- comm[i,sppInSample1] %*% t(comm[j,sppInSample2])
                    comdis[i,j] <- weighted.mean(sample.dis,sample.weights,na.rm=TRUE)
                }
                else {
                    comdis[i,j] <- mean(sample.dis,na.rm=TRUE)                
                }
            }
            else{
                comdis[i,j] <- NA
            }
        }
    }
    rownames(comdis) <- colnames(comdis) <- rownames(comm)
    as.dist(t(comdis))
}

#comdistnn
#mean distance to closest relative between taxa from two communities
comdistnt <- function(comm,dis,abundance.weighted=FALSE,exclude.conspecifics=FALSE) {
    N <- dim(comm)[1]
    comm <- decostand(comm,method="total",MARGIN=1)    
    comdisnt <- matrix(nrow=N,ncol=N)
    for (i in 1:(N-1)) {
        for (j in (i+1):N) {
            sppInSample1 <- names(comm[i, comm[i, ] > 0])
            sppInSample2 <- names(comm[j, comm[j, ] > 0])
            if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 1)) {
                sample.dis <- dis[sppInSample1, sppInSample2]
                if (exclude.conspecifics) {
                    sample.dis[sample.dis == 0] <- NA
                }
                sample1NT <- apply(sample.dis,1,min,na.rm=TRUE)
                sample2NT <- apply(sample.dis,2,min,na.rm=TRUE)
                if (abundance.weighted) {
                    sample1.weights <- comm[i,sppInSample1]
                    sample2.weights <- comm[j,sppInSample2]
                    sample.weights <- c(sample1.weights,sample2.weights)              
                    comdisnt[i,j] <- weighted.mean(c(sample1NT,sample2NT),sample.weights,na.rm=TRUE)
                }
                else
                {
                    comdisnt[i,j] <- mean(c(sample1NT,sample2NT))
                }
            }
            else{
                comdisnt[i,j] <- NA
            }
        }
    }
    rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
    as.dist(t(comdisnt))
}
