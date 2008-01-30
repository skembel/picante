`randomizeSampleKeepSppFreq` <-
function(samp)
{
	data.frame(apply(samp,2,sample),row.names=row.names(samp))
}

