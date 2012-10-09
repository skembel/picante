## draws a phylogeny with tips and nodes positioned
## by trait value, and brownian ancestral value
## node height by branch lengths
traitgram = function(
	x, phy,
	xaxt='s',
	underscore = FALSE,
	show.names = TRUE,
	show.xaxis.values = TRUE,
	method = c('ML','pic'),
	...) 
{

	method <- match.arg(method)
	Ntaxa = length(phy$tip.label)
	Ntot = Ntaxa + phy$Nnode
	phy = node.age(phy)
	ages = phy$ages[match(1:Ntot,phy$edge[,2])]
	ages[Ntaxa+1]=0
	#ages

	if (class(x) %in% c('matrix','array')) {
		xx = as.numeric(x)
		names(xx) = row.names(x)		
		} else xx = x
	
	if (!is.null(names(xx))) {
		umar = 0.1
		if (!all(names(xx) %in% phy$tip.label)) {
			print('trait and phy names do not match')
			return()
			}
		xx = xx[match(phy$tip.label,names(xx))]
		} else umar = 0.1

	lmar = 0.2
	if (xaxt=='s') if (show.xaxis.values) lmar = 1 else lmar = 0.5
	xanc <- ace(xx, phy, method=method)$ace
	xall = c(xx,xanc)
	
	a0 = ages[phy$edge[,1]]
	a1 = ages[phy$edge[,2]]
	x0 = xall[phy$edge[,1]]
	x1 = xall[phy$edge[,2]]

	tg = par(bty='n',mai=c(lmar,0.1,umar,0.1))
	if (show.names) {
		maxNameLength = max(nchar(names(xx)))
		ylim = c(min(ages),max(ages)*(1+maxNameLength/50))
		if (!underscore) names(xx) = gsub('_',' ',names(xx))
		} else ylim = range(ages)
	plot(range(c(x0,x1)),range(c(a0,a1)),
		type='n',xaxt='n',yaxt='n',
		xlab='',ylab='',bty='n',ylim=ylim,cex.axis=0.8)
	if (xaxt=='s') if (show.xaxis.values) axis(1,labels=TRUE) 
		else axis(1,labels=FALSE)
	segments(x0,a0,x1,a1)
	if (show.names) {
		text(sort(xx),max(ages),
			labels = names(xx)[order(xx)],
			adj = -0,
			srt=90)
		}
	on.exit(par(tg))
}


