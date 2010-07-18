#Evolutionary distinctiveness by: 
#a) equal splits (Redding and Mooers 2006)
#b) fair proportions (Isaac et al., 2007)
#The scale option refers to whether or not the phylogeny should be scaled to a depth of 1 or, in the case of an ultrametric tree,  scaled such that branch lengths are relative.
#If use.branch.lengths=FALSE, then all branch lengths are changed to 1.

evol.distinct <- function(tree, type=c("equal.splits", "fair.proportion"),
                            scale=FALSE, use.branch.lengths=TRUE){

type <- match.arg(type)

if(is.rooted(tree)==FALSE)
warning("A rooted phylogeny is required for meaningful output of this function", call.=FALSE)

if(scale==TRUE){
#Scale tree to have unit depth (for an ultrametric tree) or scale all branches to unit length (for an additive tree)

if(is.ultrametric(tree)==TRUE)
tree$edge.length<- tree$edge.length/(as.numeric(branching.times(tree)[1])) else 
tree$edge.length<- tree$edge.length/sum(tree$edge.length)
}

if(use.branch.lengths==FALSE)
tree$edge.length<- rep(1, length(tree$edge.length))


for(i in 1:length(tree$tip.label)){
	spp<- tree$tip.label[i]
	nodes<- .get.nodes(tree, spp)
	#get rid of root node
	nodes<- nodes[1:(length(nodes)-1)]
	
	internal.brlen<- tree$edge.length[which(tree$edge[,2] %in% nodes)]

#apportion internal branch lengths appropriately
if(length(internal.brlen)!=0){
   internal.brlen <- internal.brlen * switch(type, equal.splits = sort(rep(0.5, 
        length(internal.brlen))^c(1:length(internal.brlen))), 
        fair.proportion = {for (j in 1:length(nodes)) {
          sons <- .node.desc(tree, nodes[j])
          n.descendents <- length(sons$tips)
          if (j == 1) 
            portion <- n.descendents
          else portion <- c(n.descendents, portion)
        }; 1/portion})
}
	
	#sum internal branch lengths with the pendant edge
	ED<- sum(internal.brlen, tree$edge.length[which.edge(tree, spp)])
	
	if(i==1)
	w<- ED else
	w<- c(w, ED)
	}
results<- cbind(tree$tip.label, as.data.frame(w))
names(results)<- c("Species", "w")
return(results)
	
	}
	