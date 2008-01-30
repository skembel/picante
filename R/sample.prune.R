`sample.prune` <-
function (samp, phylo) 
{
    treeTaxa <- phylo$tip.label
    sampleTaxa <- colnames(samp)
    trimTaxa <- setdiff(treeTaxa, sampleTaxa)
    if (length(trimTaxa) > 0) drop.tip(phylo, trimTaxa) else phylo
}

