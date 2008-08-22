`.match.tree` <- function(phy, x, taxacol, traitcol, strict = FALSE) {
	# some data input error checking, all taxa in tree and x
	# no missing data values
	stopifnot(traitcol %in% names(x), taxacol %in% names(x),
		class(x) == "data.frame", class(phy) == "phylo")
	len.tips <- length(phy$tip.label)
	len.taxa <- length(x[,taxacol])
	if (any(missing <- !(phy$tip.label %in% x[, taxacol]))) {
		stop("ERROR. phylogeny tip(s): ", phy$tip.label[missing], 
			"are missing from '", substitute(x), "'")
	}
	if (strict && any(missing <- !(phy$tip.label %in% x[, taxacol])))
		stop("ERROR. '", substitute(x), "' contains taxa:", phy$tip.label[missing], 
			"not found in '", substitute(phy), "'")
	# ensure that the order of the data frame matches the tips
	# order <- match(phylo$tip.label, x[, taxcol])
	# allow taxa names as row names col names
	# deal with a vector as well
	# what to return?  
}

`df2vec` <- function(x, colID=1) {
	vec <- x[,colID]
	names(vec) <- row.names(x)
	vec
}

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

`node.age` <-
function(phy) {
	#if (phy$edge[1,1]=='-1') rootN=-1 else rootN = phy$Nnode+2
	rootN = phy$edge[1,1]

	nEdges = nrow(phy$edge)
	
	ages=rep(NA,nEdges)
	
	for (n in 1:nEdges) {
		if (phy$edge[n,1]==rootN) anc.age=0 else {
			anc.age=ages[which(phy$edge[,2]==phy$edge[n,1])]
			}
		ages[n] = anc.age + phy$edge.length[n]
		}
	phy$ages = ages
	return(phy)
	}

`sortColumns` <-
function(x) {

x[,sort(colnames(x))]

}

`sortRows` <-
function(x) {

x[sort(rownames(x)),]

}

