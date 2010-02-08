
#modified from a function by Walton Green, Gene Hunt, Scott Wing
#found at http://bricol.net/research/extinctmodel/Rfunctions/lop.R
.node.desc<- function(tr, node)
# returns vector of decendents for a particular node
{
 ee<- tr$edge
 nT<- length(tr$tip.label)
 nN<- tr$Nnode
 d.tips<- numeric()
 d.nodes<- numeric()
    
 desc<- ee[ee[,1]==node,2]
 #print (desc)
 d.tips<- append(d.tips, desc[desc<=nT])
 d.nodes<- append(d.nodes, desc[desc>nT])
 
 if (node<=nT)	d.tips<- node
 else{
 while(!all(desc<=nT))  # true if some desc are nodes not tips
  {
  desc<- ee[ee[,1] %in% desc,2]
  #print(desc)
  d.tips<- append(d.tips, desc[desc<=nT])
  d.nodes<- append(d.nodes, desc[desc>nT])  	
  }
 }
 return(list(tips=d.tips, nodes=d.nodes, node.label=node))
}