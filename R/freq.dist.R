`freq.dist` <-
function(sppFreq, metric="rank") {

#metric can be "freq" (uses difference in species frequency)
#or "rank" (default, uses difference in species frequency rank)
sppFreq <- sppFreq[sort(rownames(sppFreq)),]
if (metric=="freq") sppFreq[["rank"]] <- NULL
if (metric=="rank") sppFreq[["freq"]] <- NULL
dist(sppFreq,method="manhattan")

}

