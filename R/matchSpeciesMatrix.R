`matchSpeciesMatrix` <-
function (x, y) 
{
    mergedFrame <- merge(x, y, all.y = TRUE)
    mergedFrame[is.na(mergedFrame)] <- 0
    mergedFrame[,colnames(x)]
}

