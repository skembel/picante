`prune.sample` <-
function (samp, phylo) 
{
    treeTaxa <- phylo$tip.label
    sampleTaxa <- colnames(samp)
    trimTaxa <- setdiff(treeTaxa, sampleTaxa)
    if (length(trimTaxa) > 0) drop.tip(phylo, trimTaxa) else phylo
}

'prune.missing' <-
function(x, phylo) {
	result <- list(NULL)
    treeTaxa <- phylo$tip.label
    traitTaxa <- names(na.omit(x[phylo$tip.label]))
    trimTaxa <- setdiff(treeTaxa, traitTaxa)
    if (length(trimTaxa) > 0) 
        result$tree <- drop.tip(phylo, trimTaxa)
    else result$tree <- phylo
	result$data <- na.omit(x[phylo$tip.label])
    result
}
