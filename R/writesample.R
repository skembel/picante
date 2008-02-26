`writesample` <-
function (community, filename = "") {
    write.table(matrix2sample(community), file = filename, append = FALSE, 
        sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

