`randomizeSampleKeepSampRichness` <-
function(samp)
{
	t(data.frame(apply(samp,1,sample),row.names=colnames(samp)))
}

