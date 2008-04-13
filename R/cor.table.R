cor.table <- function (x, cor.method = "pearson", cor.type=c("standard","contrast")) 
{
    cor.type <- match.arg(cor.type)
    if (identical(cor.type,"standard")) {
        concorr <- list()
        concorr$r <- cor(x, method = cor.method)
        concorr$df <- dim(x)[1] - 2
        t <- concorr$r * sqrt(concorr$df/(1 - concorr$r^2))
        concorr$P <- dt(t, concorr$df)
        concorr
    }
    else {
    	concorr <- list()
        concorr$r <- cor(rbind(x,x*-1),method=cor.method)
        concorr$df <- length(x[,1])-1
        t <- concorr$r * sqrt(concorr$df/(1-concorr$r^2))
        concorr$P <- dt(t,concorr$df)
        concorr
    }
}
