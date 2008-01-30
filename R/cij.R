`cij` <-
function(x) {
	#Schoener index of co-occurrence
	x <- decostand(x,method="total",MARGIN=2)
	cij <- dist(t(x),method="manhattan")
	1 - (0.5 * cij)
}
