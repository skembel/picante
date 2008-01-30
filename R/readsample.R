`readsample` <-
function (filename = "") {
    x <- read.table(file = filename, header = F, sep = "\t", col.names = c("plot", 
        "abund", "id"))
    sample2matrix(x)
}

