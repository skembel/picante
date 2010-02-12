#Taxic diversity: Vane-Wright et al., 1991 and May 1990 which accounts for polytomies by counting the number of branches descending from each node that lies on the path from a spp tip to the root (not just counting the number of nodes.

tax.distinctiveness<- function(tree, type=c("Vane-Wright", "May")) {

    type <- match.arg(type)

    if(is.rooted(tree)==FALSE) {
        warning("A rooted phylogeny is required for meaningful output of this function", call.=FALSE)
    }

    n.nodes<- switch(type, "Vane-Wright"=.node.number(tree),
	"May"={edge<- tree$edge

	for(i in 1:length(tree$tip.label)){
	spp<- tree$tip.label[i]
	nodes<- .get.nodes(tree, spp)
	n.branch<- nrow(edge[edge[,1] %in% nodes,])
	if(i==1)
		n.branch.all<- n.branch else
		n.branch.all<- c(n.branch.all, n.branch)}
	n.nodes<- cbind(tree$tip.label, as.data.frame(n.branch.all))
		})

	w<- as.data.frame(1/n.nodes[,2])/sum(1/n.nodes[,2])
	results<- cbind(tree$tip.label, w)
	names(results)<- c("Species", "w")
	results
}



