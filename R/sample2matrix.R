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

