r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()

s <-  system.file("extdata", "ClusterRandomSamples.csv", package = "espatsmo") |>  read.csv()

p <- system.file("extdata", "pointsCluster.csv", package = "espatsmo") |>  read.csv()

pr <- p[s$Samples, ]

model <- ppmLGCP(points= pr, 
                 covariates = r, 
                 formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
                 dist.ar = FALSE,
                 weight.units = "km",
                 coordinates = "m")

model$base.predictions |> plot()
model$predictions |> plot()

latente <- model$predictions$mean/model$base.predictions

plot(latente)

mod.pois <- ppmSingleFit(points= pr, 
                         covariates = r, 
                         formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)",
                         as.ppmSingle = FALSE)
iml<- imFromStack(r)

pred.pois <- predict.ppm(object = mod.pois, covariates = iml, locations = iml$bio1) |> rast()

plot(pred.pois)

### Modelo cluster con corrección de sesgo
bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |>  terra::rast()

sb <- system.file("extdata", "ClusterBiasSamples.csv", package = "espatsmo") |>  read.csv()

pb <- p[sb$Samples, ] 

model.w <- ppmLGCP(points= pb, 
                   covariates = r, 
                   formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
                   dist.ar = FALSE,
                   bias.data = bias,
                   bias.correction = "weights",
                   weight.units = "km",
                   coordinates = "m")