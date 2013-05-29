#phyEstimate.R
#use phylogeny to predict trait values for new species
phyEstimate <- function(phy, trait, method="pic", ...) {
    
    #trait should be a data.frame or vector with (row)names matching phylogeny
    if (is.vector(trait)) {
        trait <- data.frame(trait)
    }
    
    trait.orig <- trait
    
    #given a tree with a novel species on it
    sppObs <- row.names(trait)
    #(novel spp. are in tree but have no trait value)
    sppUnobs <- phy$tip.label[!(phy$tip.label %in% sppObs)]
    
    res <- as.data.frame(matrix(nrow=length(sppUnobs), ncol=2, dimnames=list(sppUnobs, c("estimate","se"))))
    
    for (i in sppUnobs) {
    
        #for each novel species, prune all but measured + that species
        tree <- drop.tip(phy, subset(sppUnobs, sppUnobs != i))
    
        #root the tree at the novel species (leave root as trichotomy)
        tree <- root(tree, i, resolve.root=FALSE)
    
        #record branch length leading to novel species in rerooted tree
        edge <- Nnode(tree) - 1 + which(tree$tip.label==i)
        bl <- tree$edge.length[edge]
        
        #prune novel species and match new pruned tree <-> trait data
        tree <- drop.tip(tree, i)
        trait <- trait.orig[tree$tip.label,]
        
        #use PIC framework to estimate trait value at root node + error
        est <- ace(trait, tree, method=method, ...)
        val <- est$ace[1]
        cimax <- est$CI95[1,2]
        se <- abs(cimax-val)/1.96

        se.adj <- sqrt(bl)+se
        
        res[i,] <- data.frame(estimate=val, se=se.adj)
    
    }
    
    return(res)

}


# for discrete traits
phyEstimateDisc <- function(phy, trait, best.state=TRUE, cutoff=0.5, ...) {

    #trait should be a data.frame or vector with names matching phylogeny
    if (is.vector(trait)|is.factor(trait)) {
        trait <- data.frame(trait)
    }

    trait[,1] <- factor(trait[,1])
    trait.orig <- trait

    #given a tree with a novel taxa on it (taxa with no trait value)
    sppObs <- row.names(trait)
    sppUnobs <- phy$tip.label[!(phy$tip.label %in% sppObs)]
    trtlevels <- levels(trait[,1])
    res <- as.data.frame(matrix(nrow=length(sppUnobs), ncol=length(trtlevels), dimnames=list(sppUnobs, trtlevels)))

    #estimate support for different states for each novel taxon
    for (i in sppUnobs) {

        #for each novel species, prune all but measured + that species
        tree <- drop.tip(phy, subset(sppUnobs, sppUnobs != i))

        #root the tree at the novel species (leave root as trichotomy)
        tree <- root(tree, i, resolve.root=FALSE)

        #record branch length leading to novel species in rerooted tree
        edge <- Nnode(tree) - 1 + which(tree$tip.label==i)
        bl <- tree$edge.length[edge]
        
        #prune novel species and match new pruned tree <-> trait data
        tree <- drop.tip(tree, i)
        trait <- trait.orig[tree$tip.label,]

        #calculate value at root node and impute to novel species
        est <- ace(trait, tree, type="discrete", ...)
        val <- est$lik.anc[1,]
        
        res[i,] <- val

    }

    #estimate the best-supported state for each taxon
    if (best.state) {
        beststate <- as.data.frame(matrix(nrow=dim(res)[1], ncol=2))
        colnames(beststate) <- c("estimated.state","estimated.state.support")
        rownames(beststate) <- rownames(res)

        for (i in 1:dim(res)[1]) {       
            #if >=cutoff % taxa have same label assign a consensus taxon to node
            best <- -sort(-(res[i,]))[1]
            if (best >= cutoff) {
                beststate[i,1] <- names(best)
                beststate[i,2] <- best           
            }
            else
            {
                beststate[i,1] <- NA
                beststate[i,2] <- NA
            }
        }
    }

    #return the output
    if (best.state) {
        return(cbind(as.matrix(res),beststate))
    } else {
        return(as.matrix(res))
    }

}
