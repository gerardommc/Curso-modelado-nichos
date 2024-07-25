ZeroOneNorm <- function(x){
  require(terra)
  mi <- global(x, min, na.rm = T)[, 1]
  ma <- global(x, max, na.rm = T)[, 1]
  x.norm <- (x - mi)/(ma - mi)
}
