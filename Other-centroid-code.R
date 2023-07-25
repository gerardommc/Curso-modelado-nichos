

###########################
####Comparing centroids####
###########################

spp.test <- lapply(list.files("Real-spp/Presence/", "test", full.names = T),   
                   function(x){
                     require(rgdal)
                     x1 <- read.csv(x)
                     x1.1 <- x1
                     coordinates(x1) <- c("Longitud", "Latitud")
                     proj4string(x1) <- CRS("+init=epsg:4326")
                     x2 <- spTransform(x1, CRS = CRS("+init=epsg:4269"))
                     x1.1$x <- data.frame(coordinates(x2))[, 1]
                     x1.1$y <- data.frame(coordinates(x2))[, 2]
                     return(x1.1)}
)

ppm.preds <- list(cal.cal.pred.r, cal.mel.pred.r)
for(i in 1:2){
  spp.test[[i]]$DNC <- extract(ellips[[i]]$Suitability, spp.test[[i]][, c("x", "y")])
  spp.test[[i]]$PPM <- extract(ppm.preds[[i]], spp.test[[i]][, c("x", "y")])
}


#Correlation between surfaces and abundance
cor.test(spp.test[[1]]$Abundance, spp.test[[1]]$PPM, type = "spearman")
cor.test(spp.test[[2]]$Abundance, spp.test[[2]]$PPM, type = "spearman")

cor.test(spp.test[[1]]$Abundance, spp.test[[1]]$DNC, type = "spearman")
cor.test(spp.test[[2]]$Abundance, spp.test[[2]]$DNC, type = "spearman")

# Correlation between predicitons
cal.cal.pred.stack <- stack(cal.cal.pred.r, ellips$Cal.cal$Suitability)
names(cal.cal.pred.stack) <- c("PPM", "MVE")

cal.mel.pred.stack <- stack(cal.mel.pred.r, ellips$Cal.mel$Suitability)
names(cal.mel.pred.stack) <- c("PPM", "MVE")

pdf("Real-spp/Correlation-pred-surfaces.pdf", width = 6, height = 6)
pairs(cal.cal.pred.stack, main = "Callipepla californica")
pairs(cal.mel.pred.stack, main = "Calamospiza melanocorys")
dev.off()

##Distance between centroids
#PPM centroids

cal.cal.coef <- coefficients(cal.cal.best)
cal.mel.coef <- coefficients(cal.mel.best)

cal.cal.ppm.c <- c()
cal.mel.ppm.c <- c()

for(i in 1:3){
  cal.cal.ppm.c[i] <- -(cal.cal.coef[i+1])/(2*cal.cal.coef[i + 4])
  cal.mel.ppm.c[i] <- -(cal.mel.coef[i+1])/(2*cal.mel.coef[i + 4])
}

#Distances between centroids
cal.cal.mah <- mahalanobis(cal.cal.ppm.c, 
                           center = ellips$Cal.cal$cen.cov$centroid,
                           cov = ellips$Cal.cal$cen.cov$covariance)
cal.mel.mah <- mahalanobis(cal.mel.ppm.c, 
                           center = ellips$Cal.mel$cen.cov$centroid,
                           cov = ellips$Cal.mel$cen.cov$covariance)

#
pdf("Real-spp/Cal-cal-Results.pdf", width = 12, height = 7)
par(mfrow = c(1, 2))
plot(cal.cal.pred.r, main = "PPM")
plot(spp.buffers[[1]], add = T)
points(spp[[1]]$Longitud, spp[[1]]$Latitud,
       pch = 3, 
       col = "grey30", cex = 0.75)
plot(ellips[[1]]$Suitability, main = "MVE")
plot(spp.buffers[[1]], add = T)
points(spp.test[[1]]$x, spp.test[[1]]$y,
       cex = log10(spp.test[[1]]$Abundance)/2, pch = 20, 
       col = "grey30")
legend("topright", 
       legend = paste0(round(c(max(spp.test[[1]]$Abundance), 
                               (max(spp.test[[1]]$Abundance)-min(spp.test[[1]]$Abundance))/2, 
                               min(spp.test[[1]]$Abundance)), 1)),
       pch = 20, col = "grey30",
       bty = "n", 
       pt.cex = c(log10(max(spp.test[[1]]$Abundance))/2,
                  log10((max(spp.test[[1]]$Abundance)-min(spp.test[[1]]$Abundance))/2)/2, 
                  log10(min(spp.test[[1]]$Abundance))/2+0.1))
dev.off()

pdf("Real-spp/Cal-mel-Results.pdf", width = 12, height = 7)
par(mfrow = c(1, 2))
plot(cal.mel.pred.r, main = "PPM")
plot(spp.buffers[[2]], add = T)
points(spp[[2]]$Longitud, spp[[2]]$Latitud,
       pch = 3, 
       col = "grey30", cex = 0.75)
plot(ellips[[2]]$Suitability, main = "MVE")
plot(spp.buffers[[2]], add = T)
points(spp.test[[2]]$x, spp.test[[2]]$y,
       cex = log10(spp.test[[2]]$Abundance)/2, pch = 20, 
       col = "grey30")
legend("topright", 
       legend = paste0(round(c(max(spp.test[[2]]$Abundance), 
                               (max(spp.test[[2]]$Abundance)-min(spp.test[[2]]$Abundance))/2, 
                               min(spp.test[[2]]$Abundance)), 1)),
       pch = 20, col = "grey30",
       bty = "n", 
       pt.cex = c(log10(max(spp.test[[2]]$Abundance))/2,
                  log10((max(spp.test[[2]]$Abundance)-min(spp.test[[2]]$Abundance))/2)/2, 
                  log10(min(spp.test[[2]]$Abundance))/2+0.1))
dev.off()
