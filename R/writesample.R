`writesample` <-
function (community, filename = "") {
    write.table(matrix2sample(community), file = filename, append = FALSE, 
        sep = "\t", eol = "\n", quote = F, row.names = F, col.names = F)
}

