`.evolve.trait` <-
function(phy,
	x.root=0, #root value
	sigma = 1, #brownian motion st. dev.
	trend = 0, #brownian motion trend
	bound = c(-10,10), # bounds on trait evolution
	burst = rep(1,nrow(phy$edge)), # edge-specific rate constant
	pulse = 0, #under ACDC, probability per unit time of resetting
				#sigma to original value
	mu = rep(0,nrow(phy$edge)), #ornstein-uhlenbeck mode
	theta = 1, #ornstein-uhlenbeck parameter
	gamma = 1, #age shift factor for ACDC
	gamma.reset = 0, #probability of resetting sigma under ACDC
	trait.mode = 'Brownian',
	set.minbl = 1, #length of shortest branch
	show.plots = TRUE,
	use.color = FALSE, # plots with color
	use.bl = TRUE,
	debug = FALSE) {
	# I think this code requires the tree to be in cladewise order -- pdc
	phy <- reorder(phy, 'cladewise')
	#print(mu)
	phy$edge.length.original = phy$edge.length
	minbl = min(phy$edge.length)
	phy$edge.length = round(set.minbl*phy$edge.length/minbl)

	phy$Nterm = length(phy$tip.label)
	phy$Nedge = nrow(phy$edge)
	phy$node.traits = rep(NA,(phy$Nterm+phy$Nnode))
	if (length(grep(trait.mode,'Proportional'))==1) x.root = 1
	phy$node.traits[phy$edge[1,1]] = x.root
	phy$ages = node.age(phy)$ages
	max.ht = max(vcv.phylo(phy))
	min.trait = 0
	max.trait = 0

	if (length(grep(trait.mode,'ACDC'))==1) sg = rep(sigma,phy$Nedge)

	# Variable pu = 0 for first daughter of each node and = 1
	# for second daughter; used to assign character change to
	# one daughter in Punctuational mode
	if (length(grep(trait.mode,'Punctuational'))==1) {
		pu = rep(c(0,1),phy$Nnode)
		pu = pu[rank(phy$edge[,1],ties.method='first')]
		}
		
	edges = list(NULL)
	for (i in 1:phy$Nedge) {
		edges[[i]] = matrix(NA,phy$edge.length[i]+1,2)
		if (phy$edge[i,1]==phy$Nterm+1) {
			anc.age = 0
			anc.trait = x.root
			} else {
			anc.edge = match(phy$edge[i,1],phy$edge[,2])
			anc.age = phy$ages[match(phy$edge[i,1],phy$edge[,2])]
			anc.trait = phy$node.traits[phy$edge[i,1]]
		}
		edges[[i]][,1] = anc.age:phy$ages[i]
		edges[[i]][1,2] = anc.trait
		if (length(grep(trait.mode,'Brownian'))==1) {
				edges[[i]][-1,2] = 
					anc.trait + 
					cumsum(burst[i]*rnorm(phy$edge.length[i],
					mean = trend,
					sd = sigma))
			} else if (length(grep(trait.mode,'Bounded'))==1) {
				for (j in 2:nrow(edges[[i]])) {
					repeat {
						edges[[i]][j,2] = edges[[i]][(j-1),2] +
							rnorm(1,mean=trend,sd=sigma)
						if ((edges[[i]][j,2] > bound[1]) & 
							(edges[[i]][j,2] < bound[2])) break
					}
				}
			} else if (length(grep(trait.mode,'OU'))==1) {
				for (j in 2:nrow(edges[[i]])) {
					edges[[i]][j,2] = mu[i] + 
						theta*(edges[[i]][(j-1),2]-mu[i]) +
						rnorm(1,mean=trend,sd=sigma)
				}
			} else if (length(grep(trait.mode,'MH-OU'))==1) {
				for (j in 2:nrow(edges[[i]])) {
					edges[[i]][j,2] = mu[i] + 
						theta*(edges[[i]][(j-1),2]-mu[i]) +
						rnorm(1,mean=trend,sd=sigma)
				}
				nextbranch = which(phy$edge[,1]==phy$edge[i,2])
				mu[nextbranch] = edges[[i]][nrow(edges[[i]]),2]
			} else if (length(grep(trait.mode,'ACDC'))==1) {
				sd.edge = sg[i]
				for (j in 2:nrow(edges[[i]])) {
					#varadj = gamma^(-edges[[i]][j,1])
					edges[[i]][j,2] = edges[[i]][(j-1),2] +
						rnorm(1,
						mean=trend,
						sd=sd.edge)
					sd.edge = sd.edge/gamma
					if (runif(1) < gamma.reset) sd.edge = sigma
				}
				nextbranch = which(phy$edge[,1]==phy$edge[i,2])
				sg[nextbranch] = sd.edge
			} else if (length(grep(trait.mode,	
				'Proportional'))==1) {
				edges[[i]][-1,2] = 
					anc.trait * 
					cumprod(rlnorm(
					phy$edge.length[i],
					meanlog = trend,
					sdlog = sigma))
			} else if (length(grep(trait.mode,	
				'Speciational'))==1) {
				edges[[i]][2,2] = 
					anc.trait + rnorm(1,mean = trend,sd = sigma)
				edges[[i]][c(-1,-2),2] = edges[[i]][2,2]
			} else if (length(grep(trait.mode,	
				'Punctuational'))==1) {
				edges[[i]][2,2] = 
					anc.trait + pu[i] * 
					rnorm(1,mean = trend,sd = sigma)
				edges[[i]][c(-1,-2),2] = edges[[i]][2,2]
			} else {
				print('invalid trait evolution mode')
				break
			}
		if (min(edges[[i]][,2]) < min.trait) min.trait = min(edges[[i]][,2])
		if (max(edges[[i]][,2]) > max.trait) max.trait = max(edges[[i]][,2])
		
		phy$node.traits[phy$edge[i,2]] = edges[[i]][nrow(edges[[i]]),2]
		
		#edges = list(edges,edge)
	}

	x=matrix(phy$node.traits[1:phy$Nterm],ncol=1)
	row.names(x) = phy$tip.label
	
	if (show.plots) {
		par(mfcol=c(2,1),cex.lab=1,cex.axis=1,	
			mar=c(5.1,5.1,1.1,1.1))

		phy.coph = cophenetic.phylo(phy)
		trait.dist = as.matrix(dist(
			x,diag=TRUE,upper=TRUE))
		edgecolors <- rainbow(phy$Nedge)
		if(use.color) {
			plot.phylo(phy,show.tip.label=FALSE, edge.color = edgecolors)
		} else {
			plot.phylo(phy,show.tip.label=FALSE)
		}
		## TODO scale bar often overlaps with bottom branch
		add.scale.bar(0,1,set.minbl)	
		
		plot(c(0,max.ht),c(min.trait,max.trait),type='n',xlab='time',ylab='trait value')
		
		if(use.color) {
			for (i in 1:phy$Nedge) lines(edges[[i]], col = edgecolors[i])
		} else {
			for (i in 1:phy$Nedge) lines(edges[[i]])
		}
	
	}
	#print(phy$node.traits)
	#print(phy$edge)
	x <- as.vector(x[,1])
	names(x) <- phy$tip.label
	return(x)
}

