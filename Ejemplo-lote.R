library(espatsmo)

#Importando capas raster meiante terra
r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()

#Puntos de presencia ccompletos de especie virtual
p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()

#Muestra aleatoria representativa de la sp virtual
s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()

#Extracción de la muestra representativa
pr <- p[s$Samples, ]

plotResponses(points = pr,
              covariates = r,
              save.plot = TRUE,
              plot.pars = list(name = "Respuesta.pdf", width = 5, height = 5))

# Creando la tabla de exponentes máximos para cada una de las variables. 
#Guardar y editar a mano con base en gráficos de respeusta
exps <- data.frame(Variable = names(r), Power = NA)

#Importando base de datos de exponentes
resp <-  system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()

#identificación de conjuntos de variables ortogonales
compat <- findCompatibles(covariates = r,
                          thres = 0.7,
                          max.comb = 3)

#Creación de fórmulas
forms <- getPolyFormulas(respDF = resp, 
                         compatMat = compat)

#Ajuste de los modelos candidatos en lote
models <- ppmBatchFit(points = pr,
                      covariates = r,
                      formulas = forms,
                      parallel = FALSE,
                      top.models = 3,
                      goodness.fit = "AIC")

#Extracción del "mejor" modelo

mejor <- batchBest(models, as.ppmSingle = FALSE)

#Transformando capas raster en lista de imágenes
iml <- imFromStack(r)

pred <- spatstat.model::predict.ppm(object = mejor,
                                    covariates = iml,
                                    locations = iml$bio1) |> rast(crs = "EPSG:6372")

