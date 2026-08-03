library(espatsmo)

r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()

p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()

#Muestras con sesgo hacia zonas más muestreadas
sb <- system.file("extdata", "BiasSamples.csv", package = "espatsmo") |>  read.csv()

pb <- p[sb$Samples, ]

resp <-  system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()

compat <- findCompatibles(covariates = r,
                          thres = 0.6,
                          max.comb = 3)

forms <- getPolyFormulas(respDF = resp, 
                         compatMat = compat)

#Datos de esfuerzo de muestreo en formato raster
bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()

modelsBiasWeiRas <- ppmBatchFit(points = pb,
                                covariates = r,
                                formulas = forms,
                                bias.data = bias, #Covariable de sesgo
                                bias.correction = "weights", #Corrección de sesgo mediante modificación de pesos
                                parallel = FALSE,
                                top.models = 3)

#Datos de esfuerzo de muestreo en puntos de localidades
bias.df <- system.file("extdata", "Target-group-points.csv", package = "espatsmo") |>  read.csv()

modelsBiasWeiDF <- ppmBatchFit(points = pb,
                            covariates = r,
                            formulas = forms,
                            bias.data = bias.df,
                            bias.correction = "weights",
                            parallel = FALSE,
                            top.models = 3)

### Corrección de sesgo mediante remoción de datos de entorno

modelsBiasBacRas  <- ppmBatchFit(points = pb,
                                 covariates = r,
                                 formulas = forms,
                                 bias.data = bias, #Covariable de sesgo
                                 bias.correction = "background", #Corrección de sesgo mediante modificación de pesos
                                 parallel = FALSE,
                                 top.models = 3) 

modelsBiasBacDF  <- ppmBatchFit(points = pb,
                                covariates = r,
                                formulas = forms,
                                bias.data = bias.df, #Covariable de sesgo
                                bias.correction = "background", #Corrección de sesgo mediante modificación de pesos
                                parallel = FALSE,
                                top.models = 3) 


# Extrayendo "mejores" modelos del lote

m.sesgo1 <- batchBest(modelsBiasWeiRas, as.ppmSingle = FALSE)
m.sesgo2 <- batchBest(modelsBiasWeiDF, as.ppmSingle = FALSE)
m.sesgo3 <- batchBest(modelsBiasBacRas, as.ppmSingle = FALSE)
m.sesgo4 <- batchBest(modelsBiasBacDF, as.ppmSingle = FALSE)


