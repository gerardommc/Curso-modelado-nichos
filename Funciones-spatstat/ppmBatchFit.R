ppmBatchFit <- function(points, covariates, formulas, parallel = T, cores = 2, topModels = 10){
  require(spatstat)
  require(tidyverse)
  
  source("Funciones-spatstat/imFromStack.R")
  
  imList <- imFromStack(covariates)
  names(imList) <- names(covariates)
  
  win <- as.owin(imList[[1]])
  
  points.pp <- ppp(x = points$x, y = points$y, window = win, check = F)
  Q <- quadscheme(points.pp)
  
  if(parallel){
    require(doParallel)
    registerDoParallel(cores = cores)
    
    models <- foreach(i = seq_along(formulas)) %dopar% {
      require(spatstat)
      
      form <- as.formula(formulas[i])
      m <- ppm(Q, trend = form, 
               covariates = imList,
               na.action = na.exclude)
      
      return(m)
    }
  } else {
    models <- foreach(i = seq_along(formulas)) %do% {
      form <- as.formula(formulas[i])
      m <- ppm(Q, trend = form, 
               covariates = imList,
               na.action = na.exclude)
      
      return(m)
    }
  }
  
  perf <- sapply(models, AIC) |> sort()
  top <- perf[1:topModels]
  
  ids <- which(perf %in% top) |> sort()
  
  models <- models[ids]
  
  gc(reset = T)
  
  return(models)
}
