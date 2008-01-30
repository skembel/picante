`internal2tips` <-
function(x,int.node,return.names=FALSE) {
	# x = phy object
	# int.node = number or name of internal node
	Ntaxa = length(x$tip.label)
	Nnode = x$Nnode
	if ((Ntaxa+Nnode-1)!=nrow(x$edge)) {
		print('tree structure error')
		break
	}

	# if necessary convert int.node to a node number for an internal node
	if (mode(int.node)=='character') nodes = which(x$node.label==int.node)+Ntaxa else nodes = int.node

	tips = c()
	repeat {
		nodes = x$edge[which(x$edge[,1]%in%nodes),2]
		if (length(nodes)==0) break
		tips = c(tips,nodes)
	}
	tips = tips[tips<=Ntaxa]
	if (return.names) tips = x$tip.label[tips]
	return(tips)
}

