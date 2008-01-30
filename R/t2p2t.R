`t2p2t` <-
function(phy,species) {
	#print(phy$tip.label)
	#print(traits[,1])
	#species = data.frame(species)
	t.in.p = (species %in% phy$tip.label)
	mfp = as.character(species[!t.in.p])
	p.in.s = (phy$tip.label %in% species)
	mfs = phy$tip.label[!p.in.s]
	return(list(c('missing.from.phy:',mfp),	
		c('missing.from.species:',mfs)))
	}

