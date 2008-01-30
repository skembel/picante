`phylosignal` <-
function(x,phy,reps=999) {

  if (length(x) != length(phy$tip.label)) stop(
              "Data vector and tree contain different numbers of taxa.")
  if (!setequal(names(x), phy$tip.label)) warning(
            "Taxon names in data vector and names of tip labels do not match.")

	x <- x[phy$tip.label]
	obs.var.pic.scaled = var.pic(x,phy)
	K = Kcalc(x,phy)
	
	rnd.var.pic.scaled <- vector()
		
	#significance based on tip shuffle
	for (i in 1:reps) {
		tipsh.x <- sample(x)
		names(tipsh.x) <- names(x)
		rnd.var.pic.scaled <- c(rnd.var.pic.scaled,var.pic(tipsh.x,phy))
	}

	var.pics.scaled = c(obs.var.pic.scaled,rnd.var.pic.scaled)
	var.pics.scaled.p = rank(var.pics.scaled)[1] / (reps + 1)
	var.pics.scaled.z = (obs.var.pic.scaled - mean(var.pics.scaled))/sd(var.pics.scaled)

	data.frame(K,PIC.variance.obs=obs.var.pic.scaled,PIC.variance.random=mean(var.pics.scaled),
		PIC.variance.P=var.pics.scaled.p, PIC.variance.Z=var.pics.scaled.z)
}

