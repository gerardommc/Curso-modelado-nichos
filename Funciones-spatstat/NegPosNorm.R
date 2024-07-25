NegPosNorm <- function(x){
  ma <- global(x, max)[, 1]
  x1 <- x - ma/2
  x2 <- x1/ma
  return(x2)
}
