match.phylo.comm <- function(phy, comm) {

    if(!(is.data.frame(comm) | is.matrix(comm))) {
        stop("Community data should be a data.frame or matrix with samples in rows and taxa in columns")
    }

    res <- list()

    phytaxa <- phy$tip.label
    commtaxa <- colnames(comm)

    if(is.null(commtaxa)) {
        stop("Community data set lacks taxa (column) names, these are required to match phylogeny and community data")
    }

    if(!all(commtaxa %in% phytaxa)) {
        print("Dropping taxa from the community because they are not present in the phylogeny:")
        print(setdiff(commtaxa,phytaxa))
        comm <- comm[,intersect(commtaxa,phytaxa)]
        commtaxa <- colnames(comm)
    }
    
    if(any(!(phytaxa %in% commtaxa))) {
        print("Dropping tips from the tree because they are not present in the community data:")
        print(setdiff(phytaxa,commtaxa))
        res$phy <- prune.sample(comm, phy)
    } else {
        res$phy <- phy
    }
    
    res$comm <- comm[,res$phy$tip.label]
    return(res)

}

match.phylo.data <- function(phy, data) {

    res <- list()

    phytaxa <- phy$tip.label

    if(is.vector(data)) {
        datataxa <- names(data)
    
        if(is.null(datataxa)) {
            warning("Data set lacks taxa names, these are required to match phylogeny and data. Data are returned unsorted. Assuming that data and phy$tip.label are in the same order!")
            return(list(phy=phy,data=data))
        }
    
        if(!all(datataxa %in% phytaxa)) {
            print("Dropping taxa from the data because they are not present in the phylogeny:")
            print(setdiff(datataxa,phytaxa))
            data <- data[intersect(datataxa,phytaxa)]
            datataxa <- names(data)
        }
        
        if(any(!(phytaxa %in% datataxa))) {
            print("Dropping tips from the tree because they are not present in the data:")
            print(setdiff(phytaxa,datataxa))
            res$phy <- drop.tip(phy, setdiff(phytaxa,datataxa))
        } else {
            res$phy <- phy 
        }
        
        res$data <- data[res$phy$tip.label]
        return(res)   
        
    } else if(is.data.frame(data) | is.matrix(data)) {
        dataclass <- class(data)
        data <- as.matrix(data)

        datataxa <- rownames(data)
    
        if(is.null(datataxa)) {
            warning("Data set lacks taxa (row) names, these are required to match phylogeny and data. Data are returned unsorted. Assuming that data rows and phy$tip.label are in the same order!")
            return(list(phy=phy,data=data))
        }
    
        if(!all(datataxa %in% phytaxa)) {
            print("Dropping taxa from the data because they are not present in the phylogeny:")
            print(setdiff(datataxa,phytaxa))
            data <- data[intersect(datataxa,phytaxa),]
            datataxa <- rownames(data)
        }
        
        if(any(!(phytaxa %in% datataxa))) {
            print("Dropping tips from the tree because they are not present in the data:")
            print(setdiff(phytaxa,datataxa))
            res$phy <- drop.tip(phy, setdiff(phytaxa,datataxa))
        } else {
            res$phy <- phy
        }
        
        if (dataclass == "data.frame") {
            res$data <- as.data.frame(data[res$phy$tip.label,])
        } else {
            res$data <- data[res$phy$tip.label,]
        }
        return(res)   
             
    } else {
        stop("Data must be a vector, data.frame, or matrix")
    }
}


match.comm.dist <- function(comm, dis) {

    res <- list()

    commtaxa <- colnames(comm)

     if(is.null(commtaxa)) {
        stop("Community data set lacks taxa (column) names, these are required to match distance matrix and community data")
    }

    disclass <- class(dis)
    dis <- as.matrix(dis)

    distaxa <- rownames(dis)

    if(is.null(distaxa)) {
        warning("Distance matrix lacks taxa names, these are required to match community and distance matrix. Data are returned unsorted. Assuming that distance matrix and community data taxa columns are in the same order!")
        if (disclass == "dist") {
            return(list(comm=comm,dist=as.dist(dis))) 
        } else {
            return(list(comm=comm,dist=dis))             
        }
    }

    if(!all(distaxa %in% commtaxa)) {
        print("Dropping taxa from the distance matrix because they are not present in the community data:")
        print(setdiff(distaxa,commtaxa))
        dis <- dis[intersect(distaxa,commtaxa), intersect(distaxa,commtaxa)]
        distaxa <- rownames(dis)
    }
    
    if(any(!(commtaxa %in% distaxa))) {
        print("Dropping taxa from the community because they are not present in the distance matrix:")
        print(setdiff(commtaxa,distaxa))
        res$comm <- comm[,intersect(commtaxa,distaxa)]
    } else {
        res$comm <- comm
    }
    
    if (disclass == "dist") {
        res$dist <- as.dist(dis[colnames(comm),colnames(comm)])
    } else {
        res$dist <- dis[colnames(comm),colnames(comm)]
    }
    return(res)   

}