r <- rnorm(9, mean = 0.5, 0.5)
mig = 0.15

set.seed(123)

n.pob <- 9

K <- rpois(n.pob, 5)
image(matrix(K, 3, 3))

library(raster)

neigh <- matrix(c(0, 1, 0, 1, rep(0, 5),
                  1, 0, 1, 0, 1, rep(0, 4),
                  0, 1, 0, 0, 0, 1, rep(0, 3),
                  1, 0, 0, 0, 1, 0, 1, 0, 1,
                  0, 1, 0, 1, 0, 1, 0, 1, 0,
                  0, 0, 1, 0, 1, 0, 0, 0, 1,
                  0, 0, 0, 1, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 1, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 1, 0, 1, 0), 
                nrow = 9, ncol = 9)

image(neigh)
n.neigh <- rowSums(neigh)

N0 <- c(1, rep(0, n.pob-1))

time = 25

N <- matrix(0, nrow = n.pob, ncol = time)
N[, 1] <- N0

for(i in 2:time){
  E <- N[, i-1] * mig
  
  I <- E %*% neigh / n.neigh
  
  N[, i] <- N[, i-1 ] + r * N[, i-1] * (1- N[, i-1]/K) - E + I
}

library(animation)

saveGIF(
    for(i in 1:time){image(matrix(N[, i], 3, 3), main = i)}
)

saveGIF(
    for(i in 1:9)plot(1:25, N[i, ], type = "l", main = i)
)
