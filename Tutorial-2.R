# Simulación de presencias

library(raster); library(rgdal); library(foreach); library(spatstat)

archivos <- list.files("Datos-ejemplos/", "tif", 
                       full.names = T, 
                       recursive = F)
r <- stack(archivos)

#Definiendo un centroide para simular
centroide <- cellStats(r, mean)
r.df <- data.frame(rasterToPoints(r))
covar <- cov(r.df[, 3:5])
md <- mahalanobis(r.df[, 3:5], center = centroide, cov = covar)

md.r <- rasterFromXYZ(data.frame(r.df[, 1:2], md))
md.exp <- exp(-0.5*md.r)
plot(md.exp)

### Simulando los puntos

set.seed(182)
puntos.2 <- dismo::randomPoints(mask = md.exp,
                                n = 200,
                                prob = T)
puntos.2 <- data.frame(puntos.2)
puntos.2$x <- puntos.2$x + rnorm(200, 0, 0.05)
puntos.2$y <- puntos.2$y + rnorm(200, 0, 0.05)

### Código - favorabilidad y puntos

plot(md.exp); points(puntos.2)

# Formateo para spatstat

### Cargando las funciones

source("Funciones-spatstat/imFromStack.R")
source("Funciones-spatstat/winFromRaster.R")
source("Funciones-spatstat/plotQuantIntens.R")

#formateo-Creando todos los objetos necesarios en un solo paso
r.im <- imFromStack(r)
w <- winFromRaster(r)
puntos.2.ppp <- ppp(x = puntos.2$x,
                    y = puntos.2$y,
                    window = w,
                    check = F)
Q <- pixelquad(X = puntos.2.ppp, W = as.owin(w))

# Análisis exploratorio

### Autocorrelación
K <- envelope(puntos.2.ppp, fun = Kest, nsim = 39)
plot(K)

### Respuestas a variables

plotQuantIntens(imList = r.im,
                noCuts = 5,
                Quad = Q,
                p.pp = puntos.2.ppp,
                dir = "",
                name = "Respuestas-centroide")


### Medición de correlación entre covariables
pairs(r)

### Variables *compatibles*

#Var.1 y Var.3
# Var.2 y Var.3

#`~ Var.1 + Var.3 + I(Var.1^2) + I(Var.3^2)`
#`~ Var.2 + Var.3 + I(Var.2^2) + I(Var.3^2)`

### Ajustando los modelos

m1 <- ppm(Q = puntos.2.ppp,
          trend = ~ Var.1 + Var.3 + I(Var.1^2) + I(Var.3^2),
          covariates = r.im)
m2 <- ppm(Q = puntos.2.ppp,
          trend = ~ Var.2 + Var.3 + I(Var.2^2) + I(Var.3^2),
          covariates = r.im)

### Comparando los modelos
AIC(m1); AIC(m2)

### Analizar los efectos estimados
summary(m1)

### Diagnóstico - Residuales

par(mar = c(2,2,2,2))
diagnose.ppm(m1, main = "", cex.axis = 0.25)

### Diangnóstico - Residuales
par(mar = c(2,2,2,2))
diagnose.ppm(m2, main = "", cex.axis = 0.25)

### Diagnóstico - Ripley
K1 <- envelope(m1, fun = Kest, nsim = 39)
K2 <- envelope(m2, fun = Kest, nsim = 39)

plot(K1, cex = 0.5)

plot(K2, cex = 0.5)

### Revisando la predicción
plot(m1, se = F, main = "")

### Guardando los resultados
pred <- predict(m1)
pred.r <- raster(pred)
writeRaster(pred.r, "Predicción-m1", "GTiff",
            overwrite = T)

# Modelando los efectos espaciales

###Establecer tamaño del búfer

rr <- data.frame(r=seq(1,5,by=1))
p <- profilepl(rr, Strauss, 
               puntos.2.ppp ~ Var.1 + Var.3 + I(Var.1^2) + I(Var.3^2),
               covariates = r.im, aic=T, rbord = 0.5)

plot(p, main = "")

### Para generar un modelo de interacción

m1.int <- ppm(Q = puntos.2.ppp,
              trend = ~ Var.2 + Var.3 + I(Var.2^2) + I(Var.3^2),
              covariates = r.im,
              Strauss(p$iopt), rbord = 1) #Interacción

### Efectos estimados

coef(m1)
coef(m1.int)

## Diagnóstico, en este caso diagnose.ppm no funciona
K.int <- envelope(m1.int, Kest, nsim = 39)
plot(K.int)

### Intensidad
plot(m1.int, se = F, trend = T, cif = F)

# Proceso Cox log-Gaussiano

m1.lgcp <- kppm(puntos.2.ppp,
                trend = ~ Var.2 + Var.3 + I(Var.2^2) + I(Var.3^2),
                covariates = r.im,
                clusters = "Thomas", # Tipo de cluster, LGCP para log-Gaussian Cox Process
                statistic = "K", # K de Ripley
                method = "clik2", # Contraste con K
                model = "exp") # Modelo de varianza
### Ajustando un LGCP con `spatstat`

sum.lgcp <- summary(m1.lgcp)

# Comparando con MPP
sum.lgcp$coefs.SE.CI[, 1:4]
sum.m1$coefs.SE.CI[, 1:4]

### Predicciones 
plot(m1.lgcp, what = "intensity")

### Diagnóstico
K.lgcp <- envelope(m1.lgcp, Kest, nsim = 39)
plot(K.lgcp)