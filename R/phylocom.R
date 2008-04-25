`sample2matrix` <-
function(x) {

#turns a phylocom-format sample file into a sample X species matrix
#assumes a phylocom-format sample:
#no header line data
#column 1 = plot, column 2 = abundance, column 3 = species
#to load a phylocom-format file directly, try:
#sample2matrix(read.delim2(file="FILENAME",header=F))
	colnames(x) <- c("plot","abund","id")
    y <- tapply(x$abund, list(x$plot, x$id), sum)
    y[is.na(y)] <- 0
    y
}


`readsample` <-
function (filename = "") {
    x <- read.table(file = filename, header = FALSE, sep = "\t", col.names = c("plot", 
        "abund", "id"))
    sample2matrix(x)
}


`writesample` <-
function (community, filename = "") {
    write.table(matrix2sample(community), file = filename, append = FALSE, 
        sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


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
