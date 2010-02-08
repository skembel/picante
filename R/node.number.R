#NODE.NUMBER FUNCTION
#get the number of internal nodes for each spp

.node.number<- function(tree){
	for(i in 1:length(tree$tip.label)){
	spp<- tree$tip.label[i]
	n.node<- length(.get.nodes(tree, spp))
	if(i==1)
		n.node.all<- n.node else
		n.node.all<- c(n.node.all, n.node)
	}
n.node.all<- cbind(tree$tip.label, as.data.frame(n.node.all))
names(n.node.all)<- c("Species", "n.Nodes")
n.node.all<- n.node.all
}