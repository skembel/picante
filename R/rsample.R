
#' @title rsample
#' @description This fonction generate an community data matrix
#' @author Kerby Elpenord
#' @keywords generate, community sample
#' @param maximum       the maximum value inside an cell
#' @param tree          A phylo tree object
#' @param all_NA        if TRUE, all cells will be NA
#' @param column_with_NA if insered, the cell in the column will all be NA

rsample <-function (maximum,tree, all_NA = FALSE, column_with_NA = FALSE
){
  tree.lenth <- length(tree$tip.label)
  clump1 <- sample(0:maximum, tree.lenth, replace = TRUE)
  clump2a <- sample(0:maximum, tree.lenth, replace = TRUE)
  clump2b <- sample(0:maximum, tree.lenth, replace = TRUE)
  clump4 <- sample(0:maximum, tree.lenth, replace = TRUE)
  even <- sample(0:maximum, tree.lenth, replace = TRUE)
  random <- sample(0:maximum,tree.lenth, replace = TRUE)
  samp <- data.frame(
    clump1,
    clump2a,
    clump2b,
    clump4,
    even,
    random
  )
  samp <- t(samp)

  if (all_NA == TRUE)
  {
    samp[]<-NA
  }
  if(column_with_NA){
    samp[, column_with_NA] <- NA
  }
  colnames(samp) <- tree$tip.label
  samp
}

#' @title datasetsVerification
#' @description This fonction test an dataset for commons errors
#' @param dataset the dataset to verify
#' @details This fonction look for if it :
#' \itemize{
#' is an data.frame
#' is empty
#' contain NA
#' contain characters
#' }
#'


datasetsVerification <- function(dataset){
  isGood <- TRUE
  if (!is.data.frame(dataset)){
    print("The dataset is not an dataframe")
    isGood <- FALSE
  }
  if (all(dim(df) == c(0, 0))) {
    print("The dataset is empty")
    isGood <- FALSE
  }
  if (any(is.na(dataset))) {
    print("There are NAs in the dataset")
    isGood <- FALSE
  }
  if (any(sapply(df, function(x) any(is.character(x))))) {
    print("There are characters in the dataset")
    isGood <-FALSE
  }
  if (isGood) {
    print("The dataset seem correct")
  }
}





#' presence-absence converter
#' @description
#' Convert an data.frame to an dataset containing only 0,1 values
#' 
#'
#' @param dataset an dataframe to convert. The values in the dataset have to be
#' equal or superior to zero an to be an integer
#'
#' @return an dataset of only 0,1 values
#' @export presenceAbsenceConverter
#'
#' @examples
#' data(phylocom)
#' presenceAbsenceConverter(phylocom$sample)
presenceAbsenceConverter <- function(dataset){
  dataset[dataset > 0] = 1
  dataset
}




