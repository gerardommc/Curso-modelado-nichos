
#Paquetes
library(terra); library(foreach); library(spatstat)

#Carga y formateo para el tutorial
r <- aggregate(rast("Datos/Bioclim-CHELSA.tif"), 2) 
r <- round(scale(r), 2)

puntos <- read.csv("Datos/Puntos-analisis.csv")


fav.real <- rast("Datos/Fav-real.tif")
par(mfrow = c(1, 2))
plot(fav.real)
plot(fav.real); points(puntos, col = "red", pch = 20, cex = 0.05)

# Formateo para spatstat

### Cargando las funciones

source("Funciones-spatstat/imFromStack.R")
source("Funciones-spatstat/plotQuantIntens.R")
source("Funciones-spatstat/findCompatibles.R")
source("Funciones-spatstat/getPolyFormulas.R")
source("Funciones-spatstat/ppmBatchFit.R")

r.im <- imFromStack(r)
names(r.im) <- names(r)
w <- as.owin(r.im[[1]])
puntos.ppp <- ppp(x = puntos$x,
                  y = puntos$y,
                  window = w,
                  check = F)
Q <- pixelquad(X = puntos.ppp, W = as.owin(w))

# Análisis exploratorio

### Autocorrelación

K <- envelope(puntos.ppp, fun = Kest, nsim = 39)

plot(K)

### Respuestas a variables

plotQuantIntens(imList = r.im,
                noCuts = 5,
                Quad = Q,
                p.pp = puntos.ppp,
                dir = "",
                name = "Respuestas-centroide")

### Consideraciones para proponer modelos

### Identificación automática de covariables *compatibles*

compatibles <- findCompatibles(r, thres = 0.5, 
                               max.comb = 3)

### Obteniendo las fórmulas con `getPolyFormulas`

expon <- read.csv("Datos/Tabla-coefs.csv")
formulas  <- getPolyFormulas(respDF = expon, 
                             compatMat = compatibles)
formulas[1:5]

### Ajustando los modelos con ppmBatcchFit

modelos <- ppmBatchFit(points = puntos,
                       covariates = r, 
                       formulas = formulas[1:10],
                       parallel = F,
                       topModels = 5)

### Analizando el resultado

sapply(modelos, AIC)

summary(modelos[[4]])

### Análisis de residuales

K.modelo <- envelope(modelos[[1]], fun = "Kest", nsim = 39)

### Gráfica
par(mfrow = c(1, 2))
plot(K, main = "Datos")
plot(K.modelo, main = "Modelo")

### Métodos para residuales

  ### Gráfica de horizonte

par(mar = c(1.5,1,0,0), mfrow = c(1,1))
diagnose.ppm(modelos[[1]], cex = 0.25, outer = 5)

# Corrección de sesgo

### Definición de escenario de sesgo

sesgo <- rast(c("Datos/Target-group.tif", 
                "Datos/Distance-roads.tif"))
sesgo <- resample(sesgo, r)
plot(sesgo)

### Filtrado del entorno

source("Funciones-spatstat/maskBias.R")


### Ejemplos

r.mask1 <- maskBias(s = r[[1]], #Bio1
                    pres.areas = puntos, 
                    bias.lay = sesgo[[1]], #Target group
                    p.keep = 0.1, power = 1)

r.mask2 <- maskBias(s = r[[1]], #Bio1
                    pres.areas = puntos, 
                    bias.lay = sesgo[[1]], #Target group
                    p.keep = 0.05, power = 3, dis = 4)
r.mask3 <- maskBias(s = r[[1]], #Bio1
                    pres.areas = puntos, 
                    bias.lay = sesgo[[1]], #Target group
                    p.keep = 0.025, power = 4, dis = 4)

### Gráfica

par(mfrow = c(2, 2))
plot(r[[1]], main = "Normal")
plot(r.mask1, main = "p.keep = 0.1, power = 1")
plot(r.mask2, main = "p.keep = 0.05, power = 3")
plot(r.mask3, main = "p.keep = 0.025, power = 4")

### Histogramas

par(mfrow = c(2, 2))
hist(r[[1]], main = "Normal")
hist(r.mask1, main = "p.keep = 0.1, power = 1")
hist(r.mask2, main = "p.keep = 0.05, power = 3")
hist(r.mask3, main = "p.keep = 0.025, power = 4")


# Modelando la correlación espacial

### Modelos de interacción

### Para generar un modelo de interacción

# Establecer tamaño del búfer

rr <- data.frame(r=seq(0.1,0.5,by=0.1))
p <- profilepl(rr, AreaInter, 
               puntos.ppp ~ bio1 + bio12 + I(bio12^2) + 
                 bio18 + I(bio18^2) + 
                 I(bio18^3) + I(bio18^4),
               covariates = r.im, aic=F, rbord = 0.1)
### Para generar un modelo de interacción

par(mfrow = c(1,1))
plot(p, main = "")

### Para generar un modelo de interacción

m1.int <- ppm(Q = puntos.ppp,
              trend = ~ bio1 + bio12 + I(bio12^2) + 
                bio18 + I(bio18^2) + 
                I(bio18^3) + I(bio18^4),
              covariates = r.im,
              AreaInter(rr$r[p$iopt]), rbord = 0.1) #Interacción

### Efectos estimados

sum.int <- summary(m1.int)
sum.int$coefs.SE.CI[, 1:4]

### Efectos estimados - comparación

coef(modelos[[1]])

### Diangóstico

K.int <- envelope(m1.int, Kest, nsim = 39)

### Diangóstico

par(mfrow = c(1,2))
plot(K.modelo, main = "Modelo Poisson")
plot(K.int, main = "Modelo de interacción")

### Análisis de residuales

par(mar = c(1.5,1,0,0), mfrow = c(1,1))
diagnose.ppm(m1.int, cex = 0.25, outer = 5)


### Comparación con favorabilidad real
library(tidyverse)

modelo1 <- predict(modelos[[1]], ngrid = c(77, 116), type = "trend") |> rast()
modelo.int <- predict(m1.int, ngrid = c(77, 116), type = "trend") |> rast()

par(mfrow = c(1, 3), mar = c(1,1,1,1))
plot(fav.real, main = "real")
plot(modelo1, main = "Poisson")
plot(modelo.int,  main = "Interacción")

## Modelo con corrección de sesgo

r.s <- maskBias(s = r, 
                pres.areas = puntos, 
                bias.lay = sesgo[[1]], #Target group
                p.keep = 0.1, power = 1)

modelos.bias <- ppmBatchFit(points = puntos,
                            covariates = r.s,
                            formulas = formulas[1:20], topModels = 20)

mejor <- sapply(modelos.bias, AIC) |> which.min()

prediccion.sesgo <- predict(modelos.bias[[4]], 
                            window = w,
                            covariates = r.im, 
                            type = "trend",
                            se = F) |> rast()

par(mfrow = c(2, 2))
plot(fav.real, main = "Real")
plot(prediccion.sesgo, main = "Poisson con corrección")
plot(modelo1, main = "Poisson sin corrección")
plot(modelo.int,  main = "Interacción")
