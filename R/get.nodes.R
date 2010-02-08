#GET.NODES FUNCTION
#Get a vector of the nodes that are along the path from the root to spp x.  This may then be used to determine the number of nodes between the root and the tip (see node.number function) or the length of this path.
.get.nodes<- function(tree, spp){
	edge<- which.edge(tree, spp)
	nodes<- tree$edge[edge,1] 
	root.edge<- which(tree$edge[,1]==(length(tree$tip.label)+1))
    while(!(edge %in% root.edge)){
        edge<- which.edge(tree, tree$edge[edge,1])
        next.node<- tree$edge[edge,1]
        nodes<- c(nodes, next.node)
	}
	nodes
}