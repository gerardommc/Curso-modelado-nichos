library(espatsmo)

r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()

p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()

resp <-  system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()

compat <- findCompatibles(covariates = r,
                          thres = 0.7,
                          max.comb = 3)

forms <- getPolyFormulas(respDF = resp, 
                         compatMat = compat)

bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()

sb <- system.file("extdata", "BiasSamples.csv", package = "espatsmo") |>  read.csv()

pb <- p[sb$Samples, ]

modelsBiasRaster <- ppmBatchFit(points = pb,
                                covariates = r,
                                formulas = forms,
                                bias.data = bias,
                                bias.correction = "weights",
                                parallel = FALSE,
                                top.models = 4)

AIC(modelsBiasRaster)

v <- system.file("extdata", "ValidationSamples.csv", package = "espatsmo") |> read.csv()

pv <- p[v$Samples, ]

val.lote <- ppmValidate(model = modelsBiasRaster,
                        method = c("proc", "boyce"))

iml <- imFromStack(r)

preds <- lapply(modelsBiasRaster$models, function(x){predict.ppm(object = x, 
                                                                 covariates = iml,
                                                                 locations = iml$bio1) |> rast(crs = "EPSG:6372")})


proc <- lapply(1:4, function(i){partialROC(raster = preds[[i]],
                                           points = pv,
                                           save.plot = TRUE,
                                           plot.pars = list(name = paste0("pROC-unbias-", i, ".pdf"),
                                                            width = 5, height = 5))})

boyce <- lapply(1:4, function(i){boyceIndex(raster = preds[[i]],
                                            points = pv,
                                            save.plot = TRUE,
                                            plot.pars = list(name = paste0("Boyce-unbias-", i, ".pdf"),
                                                             width = 5, height = 5))})

desempeño <- data.frame(modelo = 1:4, 
                        AIC = AIC(modelsBiasRaster),
                        BIC = BIC(modelsBiasRaster),
                        Boyce = sapply(boyce, function(x){x$Boyce.index}),
                        pROC = sapply(proc, function(x){x$Ratio}))

write.csv(desempeño, "Metricas-evaluacion.csv", row.names = FALSE)


