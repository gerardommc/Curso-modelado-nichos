findCompatibles <- function(r, thres = 0.6, max.comb = 3){
  require(tidyverse); require(terra)
  
  rdf <- as.data.frame(r, xy = F) |> na.omit()
  
  cormat <- cor(rdf)
  
  lt <- lower.tri(cormat)
  
  cormatl <- cormat
  cormatl[!lt] <- NA
  cormatl[!(cormatl[]<thres & cormatl[]>-thres)] <- F
  cormatl[cormatl!=0] <- T
  cormatl[is.na(cormatl)] <- F
  
  n <- names(rdf)
  
  combs <- combn(n, max.comb)
  
  ref.size <- diag(max.comb) |> lower.tri() |> sum()
  
  library(foreach)
  lns <- foreach(i = 1:ncol(combs), .combine = c) %do% {
    sum(cormatl[combs[,i], combs[,i]])
  }
  
  ids <- which(lns == ref.size)
  
  df <- data.frame(t(combs[,ids]))
  
  names(df) <- paste0("Variable_", 1:nrow(combs))
  
  return(df)
}