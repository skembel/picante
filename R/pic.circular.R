### pic.R (2008-02-02)
###
###		Phylogenetically Independent Contrasts on circular data
###		an extension of the ape pic function
###
### Copyright 2008 Peter Cowan
###
### This file is *NOT* part of the R-package `ape'.

### Orginal copyright
### pic.R (2006-10-29)
###
###		Phylogenetically Independent Contrasts
###
### Copyright 2002-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

`pic.circular` <- 
function(x, phy, scaled = TRUE, var.contrasts = FALSE)
{
	if (!require(circular)) {stop("The 'circular' package is required for this function")}
	if (class(phy) != "phylo")
	  stop("object 'phy' is not of class \"phylo\"")
	if (is.null(phy$edge.length))
	  stop("your tree has no branch lengths: you may consider setting them equal
			to one, or using the function `compute.brlen'.")
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	if (nb.node != nb.tip - 1)
	  stop("'phy' is not rooted and fully dichotomous")
	if (length(x) != nb.tip)
	  stop("length of phenotypic and of phylogenetic data do not match")
	if (any(is.na(x)))
	  stop("the present method cannot (yet) be used directly with missing data: 
			you may consider removing the species with missing data from your 
			tree with the function `drop.tip'.")
	if(!is.circular(x))
	  stop("This version for circular data only")
	phy <- reorder(phy, "pruningwise")
	phenotype <- as.circular(rep(NA, nb.tip + nb.node), 
					units = 'radians', modulo = '2pi', 
					rotation = 'counter', zero = 0, type = 'angles', 
					template = 'none'
				)
	if (!is.null(names(x)) & all(names(x) %in% phy$tip.label)) {
		phenotype[1:nb.tip] <- x[phy$tip.label]
	} else {
		phenotype[1:nb.tip] <- x
		warning('the names of argument "x" and the names of the tip labels 
				did not match: the former were ignored in the analysis.')
	}

	contr <- var.con <- numeric(nb.node)

	bl=phy$edge.length

	for (i in seq(from = 1, by = 2, length.out = nb.node)) {
		j <- i + 1
		anc <- phy$edge[i, 1]
		des1 <- phy$edge[i, 2]
		des2 <- phy$edge[j, 2]
		sumbl <- bl[i] + bl[j]
		ic <- anc - nb.tip
		## get the differences between the decendant nodes
		tempcontr <- phenotype[des1] - phenotype[des2]
		abtemp <- abs(tempcontr)
		## sanity check on the difference
		if(abtemp > (2 * pi)) {
			stop("ERROR. contrast between ", substitute(phenotype[des1] ),
			 ' and ', substitute(phenotype[des2]), "is", substitute(tempcontr))
		}
		## ensure that contrasts record the short distance
		## if the diff is greater than pi recalculate 
		## then shorter side, while retaining directionality
		if(abtemp > pi) {
			if(phenotype[des1]  > phenotype[des2]) {
				tempcontr <- tempcontr - (2 * pi)
			} else if(phenotype[des1]  < phenotype[des2]) {
				tempcontr <- tempcontr + (2 * pi)
			} else {
				stop("Fatal Error", substitute(phenotype[des1] ),
				 ' and ', substitute(phenotype[des2]), " is ", substitute(tempcontr))}
		}
		contr[ic] <- tempcontr
		if (scaled) contr[ic] <- contr[ic]/sqrt(sumbl)
		if (var.contrasts) var.con[ic] <- sumbl
		
		## get vector components and weight by the branch lengths
		sin_des1 <- sin(as.numeric(phenotype[des1])) * bl[j]
		cos_des1 <- cos(as.numeric(phenotype[des1])) * bl[j]

		sin_des2 <- sin(as.numeric(phenotype[des2])) * bl[i]
		cos_des2 <- cos(as.numeric(phenotype[des2])) * bl[i]
	
		## calculate the ancestral node value
		phenotype[anc] <- as.circular(atan2(sin_des1 + sin_des2, cos_des1 + cos_des2), 
							type = "angles", units = "radians", template = "none", 
							rotation = "counter", zero = 0, modulo = "2pi"
						)
		k <- which(phy$edge[, 2] == anc)
		bl[k] <- bl[k] + bl[i]*bl[j]/sumbl
	
	}

	# TODO check that the var.contrasts = T results are as expected
	names(contr) <- 1:nb.node + nb.tip
	contr
}
