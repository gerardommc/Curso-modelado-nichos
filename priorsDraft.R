priors <- function(model, ...){
  
  UseMethod(generic = "priors", object = model)  
}

priors.ppmve <- function(model = NULL,
                         relax.precision = FALSE,
                         relax.factor = NULL){
  
  sts <- summary(model$model$samples)$statistics |> as.data.frame()
  
  rownms <- rownames(sts)
  
  cent.mean = sts$Mean[grep(pattern = "centroid.pres", x = rownms)]
  cent.prec = 1/(sts$SD[grep(pattern = "centroid.pres", x = rownms)])^2
  R = matrix(data = sts$Mean[grep(pattern = "tau.pres", x = rownms)],
             nrow = length(cent.mean),
             ncol = length(cent.mean),
             byrow = FALSE)
  beta.mean = sts$Mean[grep(pattern = "beta", x = rownms)]
  beta.prec = 1/(sts$SD[grep(pattern = "beta", x = rownms)])^2
  
  if(relax.precision){
    if(is.null(relax.factor)){stop("Please provide a precision relaxing factor > 1")}
    cent.prec = cent.prec/relax.factor
    for(i in 1:ncol(R)){R[i, i] <- R[i, i]/relax.factor}
    beta.prec = beta.prec/relax.factor
  }
  
  return(list(
    cent.mean = cent.mean,
    cent.prec = cent.prec,
    R = R,
    beta.mean = beta.mean,
    beta.prec = beta.prec
  ))  
}


