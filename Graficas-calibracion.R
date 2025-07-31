a = 18; h = 50; q = 1000; w = 10; c = 30

curve(q/h + a+ exp(-h *x/w)*(-q/h-a+c), 0, 1,
      xlab = "Tiempo (hr)", ylab = "Temperatura corporal (°C)")

x <- rnorm(100)

y <- list()

sds <- seq(0.1, 1.5, len = 5)

for(i in seq_along(sds)){
  y[[i]] <- x + rnorm(100, sd = sds[i])
}

pdf("Frecuencia.pdf", width = 5, height = 4.5)
par(mfrow = c(1,1))
density(x) |> plot(, xlab = "x", ylab = "Frecuencia", main = "")
abline(v = 0, lty = 2, col = "grey")
dev.off()

pdf("Correlacion.pdf", width = 9, height = 6.5)
par(mfrow=c(2, 3))
for(i in seq_along(y))plot(x, y[[i]], ylab = paste0("y-", i))
dev.off()

pdf("Dispersión.pdf", width = 4.5, height = 5)
plot(x, y[[i]], xlab = "x", ylab = "y")
dev.off()
