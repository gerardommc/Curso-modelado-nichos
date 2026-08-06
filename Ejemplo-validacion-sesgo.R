library(espatsmo)

r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()

p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()

s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()

pr <- p[s$Samples, ]

plotResponses(points = pr,
              covariates = r,
              save.plot = FALSE) 

resp <-  system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()

compat <- findCompatibles(covariates = r,
                          thres = 0.7,
                          max.comb = 3)

forms <- getPolyFormulas(respDF = resp, 
                         compatMat = compat)

models <- ppmBatchFit(points = pr,
                      covariates = r,
                      formulas = forms,
                      parallel = FALSE,
                      top.models = 3)

AIC(models)
BIC(models)
logLik(models)

mejor <- batchBest(models, as.ppmSingle = FALSE)

k <- envelope(mejor, fun = "Kest", nsim = 39)

plot(k)

particion <- spatialPartition(covariates = r,
                              points = pr, 
                              no.blocks = 50, 
                              part.criteria = "regular", 
                              mask.criteria = "regular",
                              seed = 432)

plot(particion, partition = c(0, 1))

val <- ppmValidate(model = mejor,
                   method = c("proc", "boyce"))

### sesgo ###

bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()

sb <- system.file("extdata", "BiasSamples.csv", package = "espatsmo") |>  read.csv()

pb <- p[sb$Samples, ]

#Corrección de sesgo con eliminación de datos de entorno
modelsBiasRaster <- ppmBatchFit(points = pb,
                                covariates = r,
                                formulas = forms,
                                bias.data = bias,
                                bias.correction = "background",
                                parallel = FALSE,
                                top.models = 3)

bias.df <- system.file("extdata", "Target-group-points.csv", package = "espatsmo") |>  read.csv()

#Corrección con modificación de pesos
modelsBiasDF <- ppmBatchFit(points = pb,
                            covariates = r,
                            formulas = forms,
                            bias.data = bias,
                            bias.correction = "weights",
                            parallel = FALSE,
                            top.models = 3)

#Modelos sin corrección del sesgo ara los datos sesgados
modelsBias <- ppmBatchFit(points = pb,
                          covariates = r,
                          formulas = forms,
                          parallel = FALSE,
                          top.models = 3)

#Ejemplp de la base de eliminación de datos de entorno
mb <- maskBias(covariates = r, bias.data = bias, points = pb)

mejor.pesos <- batchBest(modelsBiasRaster, as.ppmSingle = FALSE)
mejor.back <- batchBest(modelsBiasDF, as.ppmSingle = FALSE)
mejor.sesgo <- batchBest(modelsBias, as.ppmSingle = FALSE)

iml <- imFromStack(r)

pred <- predict.ppm(mejor, covariates = iml, locations = iml$bio1) |> rast()
pred.pesos <- predict.ppm(mejor.pesos, covariates = iml, locations = iml$bio1) |> rast()
pred.back <- predict.ppm(mejor.back, covariates = iml, locations = iml$bio1) |> rast()
pred.sesgo <- predict.ppm(mejor.sesgo, covariates = iml, locations = iml$bio1) |> rast()


plot(pred); points(pr, col = "red")
plot(pred.pesos); points(pr, col = "red")
plot(pred.back); points(pr, col = "red")
plot(pred.sesgo); points(pr, col = "red")

