mpdn<-function(sample, distance, abundance.weighted = FALSE)
{
  if(sum(colnames(sample)!=rownames(distance))>0)
  {
    sp.name <- intersect(colnames(sample), rownames(distance))
    sample <- sample[, match(sp.name, colnames(sample))]
    distance <- distance[match(sp.name, rownames(distance)), match(sp.name, rownames(distance))]
  }
  comt <- sample
  if(!abundance.weighted)
  {
    comt[comt>0] <- 1
    num <- rowSums(comt)
  }
  comt <- comt/rowSums(comt)
  comt <- as.matrix(comt)

  distance <- as.matrix(distance)
  comd <- comt %*% distance
  res <- comd * comt
  res <- rowSums(res)
  if(!abundance.weighted)
  {
    res <- res*(num/(num-1))
  }
  names(res)<-NULL
  res
}


