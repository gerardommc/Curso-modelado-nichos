##Load climate data and run PCA
library(spatstat)
library(raster)
library(foreach)

source("Funciones-spatstat/Spatstat-formatting-functions.R")

clim <- stack(list.files("Datos-centroide/", "tif", full.names = T))

## Load spp data and format
calcal <-  read.csv("Datos-centroide/cal_cal_test.csv")
calmel <- read.csv("Datos-centroide/cal_mel_test.csv")

spp <- list(calcal, calmel)

spp.buffers <- lapply(spp, function(x){
   require(rgdal); require(rgeos)
   coordinates(x) <- c("Longitud", "Latitud")
   buf <- rgeos::gBuffer(x, width = 5, byid = F)
   return(buf)
})

## Cropping climate
clim.spp <- lapply(spp.buffers, function(x1){
   m <- mask(x = clim, mask = x1)
   m.df <- na.omit(rasterToPoints(m))
   m <- rasterFromXYZ(m.df)
   return(m)
})

##Working windows
win <- list(
   winFromRaster(raster(clim.spp[[1]][[1]])),
   winFromRaster(raster(clim.spp[[2]][[1]]))
) 

## transforming species points to a planar point pattern
spp.ppp <- lapply(1:2, function(x){
      ppp(spp[[x]][, 'Longitud'], spp[[x]][, 'Latitud'], window = win[[x]], check = F)
})

## Transforming species' layers to spatstat images
bio.im.list <- lapply(clim.spp, imFromStack)


##Responses to bioclimatic variables
###Model formulas from compatibe and suitable variables

cal.cal.formula <- "~ bio8 + bio11 + bio12 + I(bio8^2) + I(bio11^2) + I(bio12^2)"

cal.mel.formula <- "~ bio5 + bio7 + bio16 + I(bio5^2) + I(bio7^2) + I(bio16^2)"

## Ajuste de modelo
cal.cal.ppm <- ppm(spp.ppp[[1]],
                trend = formula(cal.cal.formula),
                covariates = bio.im.list[[1]])

cal.mel.ppm <- ppm(spp.ppp[[2]],
       trend = formula(cal.mel.formula),
       covariates = bio.im.list[[2]])

#Diagnóstico de los modelos

summary(cal.cal.ppm)
summary(cal.mel.ppm)

diagnose.ppm(cal.cal.ppm)

diagnose.ppm(cal.mel.ppm)

cal.cal.pred <- predict(cal.cal.ppm, type = "intensity", dimyx = c(168, 116))
cal.mel.pred <- predict(cal.mel.ppm, type = "intensity", dimyx = c(194, 178))

cal.cal.pred.r <- raster(cal.cal.pred)
cal.mel.pred.r <- raster(cal.mel.pred)

plot(cal.cal.pred)
plot(cal.mel.pred)

########################
####Fitting ellipses####
########################
spp.vars <- list(c(8, 11, 12), c(5, 7, 16))

clim.spp.points <- lapply(clim.spp, function(x){data.frame(rasterToPoints(x))})

ellips <- foreach(i = seq_along(spp)) %do% {
   data <- data.frame(extract(clim.spp[[i]], spp[[i]]))[, paste0("bio", spp.vars[[i]])]
   data <- na.omit(data)
   cent <- colMeans(data)
   cov <- cov(data)
   dist <- mahalanobis(clim.spp.points[[i]][, names(data)], center = cent, cov = cov)
   dist.r <- rasterFromXYZ(data.frame(clim.spp.points[[i]][, c("x", "y")], dist))
   suit <- exp(-0.5 * dist.r)
   return(list(cen.cov = cent, 
               Suitability = suit,
               Distance = dist.r))
}

names(ellips) <- c("Cal.cal", "Cal.mel")

#Revisando centroides

coef.calcal <- coef(cal.cal.ppm)
coef.calmel <- coef(cal.mel.ppm)

cent.calcal <- c()
cent.calmel <- c()

for(i in 1:3){
 cent.calcal[i] <- -coef.calcal[i+1]/(2*coef.calcal[i + 4])
 cent.calmel[i] <- -coef.calmel[i+1]/(2*coef.calmel[i + 4])
}

cent.calcal
ellips$Cal.cal$cen.cov

cent.calmel
ellips$Cal.mel$cen.cov

par(mfrow  = c(1, 2))
plot(cal.cal.pred.r, main = "PPM")
plot(ellips$Cal.cal$Suitability, main = "Elipsoide")

plot(cal.mel.pred.r, main = "PPM")
plot(ellips$Cal.mel$Suitability, main = "Elipsoide")

