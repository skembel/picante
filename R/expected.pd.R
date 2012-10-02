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

# utility function - modified from ape's extract.clade function
.extract.clade.noreord <- function (phy, node, root.edge = 0, interactive = FALSE) 
{
  Ntip <- length(phy$tip.label)
  ROOT <- Ntip + 1
  Nedge <- dim(phy$edge)[1]
  wbl <- !is.null(phy$edge.length)
  if (interactive) 
    node <- identify(phy)$nodes
  else {
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
  }
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