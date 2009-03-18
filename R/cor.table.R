cor.table <- function (x, cor.method = c("pearson","spearman"),
    cor.type = c("standard", "contrast")) 
{
    cor.method <- match.arg(cor.method)
    cor.type <- match.arg(cor.type)
    if (identical(cor.type, "standard")) {
        concorr <- list()
        concorr$r <- cor(x, method = cor.method)
        concorr$df <- dim(x)[1] - 2

    }
    else {
        concorr <- list()
        concorr$r <- cor(rbind(x, x * -1), method = cor.method)
        concorr$df <- length(x[, 1]) - 1
    }
    t <- concorr$r * sqrt(concorr$df/(1 - concorr$r^2))
    concorr$P <- 2*pt(t, concorr$df)
    concorr$P[concorr$r>0]<-2*pt(t[concorr$r>0], concorr$df,lower.tail=FALSE)
    concorr
}
