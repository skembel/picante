`writetraits` <-
function (trt, file = "", bin = NULL, sigd = 3) 
{
    head = matrix(c("3"), ncol = length(names(trt)), nrow = 1)
    if (!is.null(bin)) 
        head[bin] = 0
    write.table(data.frame("type", head), file = file, sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(matrix(c("name", colnames(trt)), nrow = 1), 
        file = file, sep = "\t", append = TRUE, quote = FALSE, 
        row.names = FALSE, col.names = FALSE)
    write.table(signif(trt, sigd), sep = "\t", file = file, append = TRUE, quote = FALSE, 
        row.names = TRUE, col.names = FALSE)
}

