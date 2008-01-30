`sppFreq` <-
function(x) {

freqVec <- apply(decostand(x,"pa"),2,sum)/nrow(x)
rankVec <- rank(1-freqVec,ties.meth="min")
freqVec <- data.frame("freq"=freqVec,"rank"=rankVec)
freqVec[sort(rownames(freqVec)),]

}

