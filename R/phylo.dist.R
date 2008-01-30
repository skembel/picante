`phylo.dist` <-
function (phylo) 
{
    phylo.d <- cophenetic(phylo)
    phylo.d <- phylo.d[sort(rownames(phylo.d)), sort(colnames(phylo.d))]
    as.dist(phylo.d)
}

