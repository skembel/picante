`readsample` <-
function (filename = "") {
    x <- read.table(file = filename, header = FALSE, sep = "\t", col.names = c("plot", 
        "abund", "id"))
    sample2matrix(x)
}

