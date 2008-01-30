'pruneMissing' <-
function(x,tree) {
	result <- list(NULL)
    treeTaxa <- tree$tip.label
    traitTaxa <- names(na.omit(x[tree$tip.label]))
    trimTaxa <- setdiff(treeTaxa, traitTaxa)
    if (length(trimTaxa) > 0) 
        result$tree <- drop.tip(tree, trimTaxa)
    else result$tree <- tree
	result$data <- na.omit(x[tree$tip.label])
    result
}
