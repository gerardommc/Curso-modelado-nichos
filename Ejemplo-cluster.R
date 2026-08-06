library(espatsmo)

r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()

s <-  system.file("extdata", "ClusterRandomSamples.csv", package = "espatsmo") |>  read.csv()

p <- system.file("extdata", "pointsCluster.csv", package = "espatsmo") |>  read.csv()

pr <- p[s$Samples, ]

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
  