maskBias <- function(s, pres.areas, bias.lay, p.keep = 0.1, power = 1, dis = NA){
  require(terra)
  
  ZeroOneNorm <- function(x){
    mi <- global(x, min, na.rm = T)[, 1]
    ma <- global(x, max, na.rm = T)[, 1]
    x.norm <- (x - mi)/(ma - mi)
  }
  
  if(!is.na(dis)){s <- disagg(s, fact = dis)}
  
  s.r <- s
  bias.lay <- ZeroOneNorm(bias.lay^power)
  bias.df <- as.data.frame(bias.lay, xy = T)
  
  qs <- quantile(bias.df[, 3], probs = seq(0, 0.9, len = 9), na.rm = T)
  
  qs.r <- lapply(qs, function(x){
    q1 <- bias.lay>x
    q1 <- classify(q1, rcl = matrix(c(1,1, 0, NA), byrow = T, ncol = 2))
    b1 <- mask(bias.lay, q1)
    return(b1)
  })
  
  bg <- lapply(qs.r, function(x){ #Aquí está el probias.layema
    sam <- sample(1:nrow(bias.df), size = ceiling(nrow(bias.df)*p.keep), prob = bias.df[, 3], replace = F)
    pt <- bias.df[sam, c("x", "y")]
    pt.r <- rasterize(as.matrix(pt), bias.lay)
    return(pt.r)})
  
  bg.m <- merge(bg[[1]], bg[[2]], bg[[3]], bg[[5]], bg[[6]], bg[[7]], bg[[8]], bg[[9]])
  
  pres <- rasterize(as.matrix(pres.areas), bias.lay)
  
  pres <-  classify(pres, rcl = matrix(c(1, Inf, 1, 0, 0, NA), byrow = T, ncol = 3))
  
  m.samp <- mask(s.r, bg.m, inverse = F)
  
  p.samp <- mask(s.r, pres)
  
  s.bias <- merge(m.samp, p.samp)
  
  return(s.bias)
}

