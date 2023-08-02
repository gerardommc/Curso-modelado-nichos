r <- 0.1; mig = 0.05; M <- 5

set.seed(123)

n.pob <- 9

K <- round(rnorm(n.pob, 30, 12))
image(matrix(K, 3, 3))

neigh <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0,
                  1, 0, 1, 0, 1, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 1, 0, 0, 0,
                  1, 0, 0, 0, 1, 0, 1, 0, 0,
                  0, 1, 0, 1, 0, 1, 0, 1, 0,
                  0, 0, 1, 0, 1, 0, 0, 0, 1,
                  0, 0, 0, 1, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 1, 0, 1, 0, 1,
                  0, 0, 0, 0, 0, 1, 0, 1, 0), 
                nrow = 9, ncol = 9)

image(neigh)
n.neigh <- rowSums(neigh)

N0 <- round(rnorm(9, mean = 5, sd = 3))

time = 250

N <- matrix(0, nrow = n.pob, ncol = time)
N[, 1] <- N0

dt = 0.1

for(i in 2:time){
  E <- N[, i-1] * mig
  
  I <- E %*% neigh / n.neigh
  
  N[, i] <- N[, i-1 ] + (r * N[, i-1] * (1- N[, i-1]/K) * (N[, i-1] - M) - E + I) * dt
}

library(animation)

saveGIF(
    for(i in seq(1, 250, by = 10)){image(matrix(N[, i], 3, 3), main = i)}
)

saveGIF(
    for(i in 1:9)plot(1:time, N[i, ], type = "l", main = i)
)
