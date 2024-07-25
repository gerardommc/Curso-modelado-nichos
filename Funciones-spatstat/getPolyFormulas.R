getPolyFormulas <- function(respDF, compatMat){
  require(foreach)
  
  formulas <-  foreach(i = 1:nrow(compatMat), .combine = c) %do% {
    v <- compatMat[i, ]
    rd <- respDF[which(respDF$Variable %in% v), ]
    
    exponents <- foreach(j = 1:nrow(rd)) %do% {
      p <- 1:rd$Power[j]
    }
    
    var.exponents <- foreach(k = seq_along(exponents), .combine = c) %do% {
      ifelse(exponents[[k]] == 1, v[k], paste0("I(", v[k], "^", exponents[[k]], ")"))
    }
    
    f1 <- paste0("~",var.exponents[1])
    for(ii in 2:length(var.exponents)){
      f1 <- paste(f1, var.exponents[ii], sep = " + ")
    }
    return(f1)
  }
  
  return(formulas) 
}


