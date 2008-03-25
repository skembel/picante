`internal2tips` <-
function(phy,int.node,return.names=FALSE) {
	# phy = phy object
	# int.node = number or name of internal node
	Ntaxa = length(phy$tip.label)
	Nnode = phy$Nnode
	if ((Ntaxa+Nnode-1)!=nrow(phy$edge)) {
		print('tree structure error')
		break
	}

	# if necessary convert int.node to a node number for an internal node
	if (mode(int.node)=='character') nodes = which(phy$node.label==int.node)+Ntaxa else nodes = int.node

	tips = c()
	repeat {
		nodes = phy$edge[which(phy$edge[,1]%in%nodes),2]
		if (length(nodes)==0) break
		tips = c(tips,nodes)
	}
	tips = tips[tips<=Ntaxa]
	if (return.names) tips = phy$tip.label[tips]
	return(tips)
}

