#' @name cor.table
#' @title  Table of correlations and P-values
#' @description
#' Table of correlations with associated P-values and df,
#'  can be used with regular or independent contrast data
#' @usage cor.table(x, cor.method = c("pearson","spearman"),
#' cor.type=c("standard","contrast")) 
#' @param x Data frame of data points or contrasts at nodes
#' @param cor.method  Correlation method
#' @param cor.type Are data standard or independent contrast values?
#' @keywords univar 
#' @references Garland, T., Jr., P. H. Harvey, and A. R. Ives. 1992.
#'  Procedures for the analysis of comparative data using phylogenetically
#'   independent contrasts. Systematic Biology 41:18-32.
#' Table of correlations and P-values
#' Table of correlations with associated P-values and df, can be used with
#' regular or independent contrast data
#' @examples
#' example_1 <- data.frame(
#' X = c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968),
#' Y = c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259))
#' cor.table(x          = example_1, 
#'          cor.method = "pearson",
#'          cor.type   = "contrast")
#' @export cor.table
cor.table <- function(
    x, cor.method = c("pearson", "spearman"),
    cor.type = c("standard", "contrast")) {
  cor.method <- match.arg(cor.method)
  cor.type <- match.arg(cor.type)
  if (identical(cor.type, "standard")) {
    concorr <- list()
    concorr$r <- cor(x, method = cor.method)
    concorr$df <- dim(x)[1] - 2
  } else {
    concorr <- list()
    concorr$r <- cor(rbind(x, x * -1), method = cor.method)
    concorr$df <- length(x[, 1]) - 1
  }
  t <- concorr$r * sqrt(concorr$df / (1 - concorr$r^2))
  concorr$P <- 2 * pt(t, concorr$df)
  concorr$P[concorr$r > 0] <- 2 * pt(t[concorr$r > 0],
    concorr$df,
    lower.tail = FALSE
  )
  concorr
}
