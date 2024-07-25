regRaster <- function(s){
  av <- global(s, mean, na.rm = T)[, 1]
  desv <- global(s, sd, na.rm = T)[, 1]
  s1 <- (s-av)/desv
  return(s1)
}
