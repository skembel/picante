`randomizeSpeciesMatrix` <-
function (x, keepSppFreq = TRUE) 
{
    if (keepSppFreq) sppFreq <- apply(decostand(x, "pa"), 2, sum)/ncol(x)
	    else sppFreq <- apply(decostand(x, "pa"), 2, max)/ncol(x)
    siteRichness <- apply(decostand(x, "pa"), 1, sum)
    sampleList <- vector()
    sppList <- vector()
    for (siteNum in 1:length(siteRichness)) {
        sampleList <- c(sampleList, rep(names(siteRichness)[siteNum], 
            siteRichness[siteNum]))
        sppList <- c(sppList, sample(names(sppFreq), siteRichness[siteNum], 
            replace = FALSE, prob = sppFreq))
    }
    shuffledList <- data.frame(sample = sampleList, species = sppList, 
        p = rep(1, length(sppList)))
    shuffledMatrix <- tapply(shuffledList$p, list(shuffledList$sample, 
        shuffledList$species), sum)
    shuffledMatrix[is.na(shuffledMatrix)] <- 0
    shuffledMatrix[rownames(x), ]
}

