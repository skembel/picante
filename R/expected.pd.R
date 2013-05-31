# calculate expected PD for all subsets of a phylogeny
expected.pd <- function(phy) {
  ead <- ead(phy)
  epd <- vector(mode="numeric", length=Ntip(phy))
  for(n in 1:Ntip(phy)){
    epd[n] <- sum(ead$edge.length * (1 - dbinom(0, ead$num.children, n/Ntip(phy))))
  }
  return(data.frame(n=1:Ntip(phy),expected.pd=epd))
}

# calculate the edge abundance distribution (EAD)
ead <- function(phy) {
  phy <- reorder(phy)
  Nedges <- length(phy$edge.length)
  edgelen <- vector(mode = "numeric", length=Nedges)
  children <- vector(mode = "numeric", length=Nedges)
  for (i in 1:Nedges) {
    edgelen[i] <- phy$edge.length[i]
    descnode <- phy$edge[i,2]
    if (descnode <= Ntip(phy)) {
      children[i] <- 1
    }
    else
    {
      children[i] <- Ntip(.extract.clade.noreord(phy, descnode))
    }
  }
  rawdata <- data.frame(num.children=children, edge.length=edgelen)
  byclass <- aggregate(rawdata$edge.length, by=list(num.children=rawdata$num.children), sum)
  colnames(byclass)[2] <- "edge.length"
  return(byclass)
}

#expected variance of PD
variance.pd <- function(phy, upper.bound=TRUE) {
  Nedges <- length(phy$edge.length)
  if (upper.bound) {
    ead <- .edge.children.var.ub(phy)
    epd <- .binomial.pd(ead,Ntip(phy))
    varpd <- .var.pd.ub(ead,Ntip(phy))  
  } else {
    ead <- .edge.children.var.ub(phy)
    eadvar <- .edge.children.var(phy)
    epd <- .binomial.pd(ead,Ntip(phy))
    varpd <- .var.pd(eadvar,Ntip(phy))
  }
  return(data.frame(expected.pd=epd, variance.pd=varpd))
}

# binomial PD
.binomial.pd <-function(ead,Ntip){
  epd=vector(mode="numeric",length=Ntip)
  for(n in 1:Ntip){
    epd[n]<-sum(ead$byclass$edge.length*(1-dbinom(0,ead$byclass$num.children,n/Ntip)))
  }
  return(epd)
}

#expected variance in pd from R(k),R(k,l)
.var.pd <- function(test, Ntip) {
  varpd=vector(mode="numeric",length=Ntip)
  for(n in 1:Ntip) {
  	varpd[n]=sum(test$byclass.sq$sum.sq.edge.length*((1-n/Ntip)^test$byclass.sq$num.children)*(1-(1-n/Ntip)^test$byclass.sq$num.children))+2*sum(test$byclass.prod$sum.prods.edge.length*((1-n/Ntip)^test$byclass.prod$num.children)*(1-(1-n/Ntip)^test$byclass.prod$num.desc.children))  
	}
	return(varpd)
}

# upper bound on variance in PD---quicker by far than var.pd
.var.pd.ub<-function(test,Ntip) {
  varpd <- vector(mode="numeric",length=Ntip)
  for(n in 1:Ntip) {
	  varpd[n]=sum((test$byclass.sq$sum.sq.edge.length+2*test$byclass.prod$edge.length.times.subclade)*((1-n/Ntip)^test$byclass.sq$num.children)*(1-(1-n/Ntip)^test$byclass.sq$num.children))
  }
  return(varpd)
}

# edge children variance
.edge.children.var <- function(phy) {
  Nedges <- length(phy$edge.length)
  edgelen <- vector(mode = "numeric", length=Nedges)
  children <- vector(mode = "numeric", length=Nedges)
  rawdata<-data.frame()
  for (i in 1:Nedges) {  
    edgelen[i] <- phy$edge.length[i]
    descnode <- phy$edge[i,2]
    if (descnode <= Ntip(phy)) {
      children[i] <- 1
    }
    else
    {
      desc.phy<-.extract.clade.noreord(phy, descnode)
      children[i] <- Ntip(desc.phy)
      desc.Nedges <- length(desc.phy$edge.length)
      desc.edgelen <- vector(mode = "numeric", length=desc.Nedges)
      desc.children <- vector(mode = "numeric", length=desc.Nedges)
      for (j in 1:desc.Nedges) {    
        desc.edgelen[j] <- desc.phy$edge.length[j]
        desc.descnode <- desc.phy$edge[j,2]
        if (desc.descnode <= Ntip(desc.phy)) {
          desc.children[j] <- 1
        }
        else
        {
          desc.children[j] <- Ntip(.extract.clade.noreord(desc.phy, desc.descnode))
        }
      } 
      rawdata <- rbind(rawdata,data.frame(num.children=matrix(children[i],nrow=length(desc.children), ncol=1), num.desc.children=desc.children, prod.edge.length=edgelen[i]*desc.edgelen)) 
    }
  }
  byclass.prod <- aggregate(rawdata$prod.edge.length, by=list(num.children=rawdata$num.children,num.desc.children=rawdata$num.desc.children), sum)        
  rawdata2<- data.frame(num.children=children, node.edge.length=edgelen)  
  byclass.sq <-   aggregate(rawdata2$node.edge.length^2, by=list(num.children=rawdata2$num.children), sum)     
  colnames(byclass.prod)[3] <- "sum.prods.edge.length"
  colnames(byclass.sq)[2] <- "sum.sq.edge.length"
  return(list(rawdata=rawdata, byclass.prod=byclass.prod,byclass.sq=byclass.sq))
}

# calculate components needed for upper bound on variance
.edge.children.var.ub <- function(phy) {
  phy <- reorder(phy)
  Nedges <- length(phy$edge.length)
  edgelen <- vector(mode = "numeric", length=Nedges)
  children <- vector(mode = "numeric", length=Nedges)
  product<-vector(mode = "numeric", length=Nedges)
  for (i in 1:Nedges) {
    edgelen[i] <- phy$edge.length[i]
    descnode <- phy$edge[i,2]
    if (descnode <= Ntip(phy)) {
      children[i] <- 1
    }
    else
    {
      desc.tree<-.extract.clade.noreord(phy, descnode)  
      children[i] <- Ntip(desc.tree)
      product[i]<-edgelen[i]*(sum(desc.tree$edge.length))
    }
  }
  rawdata <- data.frame(num.children=children, edge.length=edgelen,product=product)
  byclass <- aggregate(rawdata$edge.length, by=list(num.children=rawdata$num.children), sum)
  byclass.sq <-   aggregate(rawdata$edge.length^2, by=list(num.children=rawdata$num.children), sum)
  byclass.prod<-   aggregate(rawdata$product, by=list(num.children=rawdata$num.children), sum)
  colnames(byclass)[2] <- "edge.length"
  colnames(byclass.sq)[2] <- "sum.sq.edge.length"
  colnames(byclass.prod)[2] <- "edge.length.times.subclade"
  return(list(rawdata=rawdata, byclass=byclass,byclass.sq=byclass.sq,byclass.prod=byclass.prod))
}

# utility function - modified from ape's extract.clade function
.extract.clade.noreord <- function (phy, node, root.edge = 0) 
{
  Ntip <- length(phy$tip.label)
  ROOT <- Ntip + 1
  Nedge <- dim(phy$edge)[1]
  wbl <- !is.null(phy$edge.length)
	if (length(node) > 1) {
	  node <- node[1]
	  warning("only the first value of 'node' has been considered")
	}
	if (is.character(node)) {
	  if (is.null(phy$node.label)) 
		stop("the tree has no node labels")
	  node <- which(phy$node.label %in% node) + Ntip
	}
	if (node <= Ntip) 
	  stop("node number must be greater than the number of tips")
  if (node == ROOT) 
    return(phy)
  # phy <- reorder(phy)
  root.node <- which(phy$edge[, 2] == node)
  start <- root.node + 1
  anc <- phy$edge[root.node, 1]
  next.anc <- which(phy$edge[-(1:start), 1] <= anc)
  keep <- if (length(next.anc)) 
    start + 0:(next.anc[1] - 1)
  else start:Nedge
  if (root.edge) {
    NewRootEdge <- phy$edge.length[root.node]
    root.edge <- root.edge - 1
    while (root.edge) {
      if (anc == ROOT) 
        break
      i <- which(phy$edge[, 2] == anc)
      NewRootEdge <- NewRootEdge + phy$edge.length[i]
      root.edge <- root.edge - 1
      anc <- phy$edge[i, 1]
    }
    if (root.edge && !is.null(phy$root.edge)) 
      NewRootEdge <- NewRootEdge + phy$root.edge
    phy$root.edge <- NewRootEdge
  }
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  TIPS <- phy$edge[, 2] <= Ntip
  tip <- phy$edge[TIPS, 2]
  phy$tip.label <- phy$tip.label[sort(tip)]
  phy$edge[TIPS, 2] <- order(tip)
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[sort(unique(phy$edge[,1])) - Ntip]
  Ntip <- length(phy$tip.label)
  phy$Nnode <- dim(phy$edge)[1] - Ntip + 1L
  newNb <- integer(Ntip + phy$Nnode)
  newNb[node] <- Ntip + 1L
  sndcol <- phy$edge[, 2] > Ntip
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (Ntip + 
    2):(Ntip + phy$Nnode)
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  phy
}