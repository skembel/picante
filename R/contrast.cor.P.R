`contrast.cor.P` <-
function(r,df) {
	t <- r * sqrt(df/(1-r^2))
	dt(t,df)
}

