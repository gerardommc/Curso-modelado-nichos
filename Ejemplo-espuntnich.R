library(espuntnich)

# Cargando lo datos ambientales
r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |> terra::rast() |> scale()

#Los puntos de presencia
p <- system.file("extdata", "points.csv", package = "espatsmo") |> read.csv()

s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()

#Extracción de la muestra representativa
pr <- p[s$Samples, ]
## Llamando a el modelo donde la intensidad de puntos es una función de 
# la distancia mahalanobis al centroide del nicho.

m <- ppmve(points = pr,
           covariates = r,
           covariate.names = names(r), #Nombres de las covariables
           CovMat = "local", #Especificar si la matriz de covarianza se parametrizará mediante las condiciones en los sitios de presencia (local)
           Distance = "mahalanobis", #Tipo de distancia, la alternativa es euclidean
           no.bkgd = 5000, # Número de puntos de entorno
           niter = 10000,
           nthin = 9,
           nburnin = 1000,
           chains = 2,
           parallel = TRUE) # Poner parallel = FALSE si su computadora sólo tiene 1 núcleo

m$model$samples |> coda::traceplot()

summary(m)

preds <- predict(object = m, newdata = r, probs = c(0.25, 0.5, 0.975))

plot(preds)
plot(preds$Prob.0.5); points(p, col = "red")

resumen <- summary(m)

prec <- matrix(data = c(resumen$statistics[10:25, 1]), nrow = 4, ncol = 4, byrow = FALSE)

pr <- mahalLocallocalPriors(cent.mean = resumen$statistics[2:5, 1],
                              cent.prec = 1/(resumen$statistics[2:5, 2])^2,
                              R = prec,
                              beta.mean = resumen$statistics[1, 1],
                              beta.prec = 1/(resumen$statistics[1, 2])^2)

m.loc <-ppmve(points = p,
              covariates = r,
              covariate.names = names(r),
              CovMat = "locallocal", #Para obtener muestras de la matriz de covarianza mediante la intensidad de puntos
              Distance = "mahalanobis",
              priors = pr,
              no.bkgd = 5000,
              niter = 10000,
              nthin = 9,
              nburnin = 1000,
              chains = 2, 
              parallel = TRUE)


m.loc$model$samples |> coda::traceplot()
m.loc$model$samples |> coda::geweke.diag() #Diagnóstico de convergencia

summary(m.loc)

preds.loc <- predict(object = m.loc, newdata = r, probs = c(0.25, 0.5, 0.975))

plot(preds.loc)
plot(preds)

