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

