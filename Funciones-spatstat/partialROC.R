#This function calculates the Partial ROC statistic (Peterson et al. 2008) using as inputs:
# raster = A raster layer object in format terra:SpatRaster
# points = A two-column (accepts three columns) data frame where the first columns have to be the x-coordinates and the second are the y-coordinates
# p.points = the proportion of validation points used in each iteration
# iterations = the number of times the random sampling is repeated
# pdf.name = Name of the file to save a graph of all the curves
# buf = a radius around testing presence points which will eliminate all areas further away (redundant if omission > 0)
# omission = quantile of suitability values which will be used to exclude presence training points (my personal interpretation of Peterson et al. 2008).
### Ideally when running the test the result has to be stored in an R object which will contain the different areas and the ratios. For significance simply calculate the proportion of Ratios < 1. 
### Issues: With small sample sizes, ratios can occaionally be > 2 because the randomised random expectation area could be << 1

partialROC <- function(raster, points, p.points = 0.5, iterations = 39,
                       buf = NULL, omission = 0.01, save.pdf = T, pdf.name = "PartialROC"){
  require(terra)
  
  ZeroOneNorm <- function(x){
    require(terra)
    mi <- global(x, min, na.rm = T)[, 1]
    ma <- global(x, max, na.rm = T)[, 1]
    x.norm <- (x - mi)/(ma - mi)
  }

  if(!is.null(buf)){
    p <- vect(as.matrix(points[, -3]))
    bu <- buffer(p, width = buf)
    
    bur <- rasterize(bu, raster)
    
    raster <- mask(raster, bur)
  }
  
  raster <- log(raster + 0.1)
  
  if(omission > 0){
    vals <- terra::extract(raster, points[, 1:2])[,2]
    q <- quantile(vals, omission, na.rm = T)
    om.r <- raster > q
    om.r <- classify(om.r, rcl = matrix(c(-Inf, 0, NA), ncol = 3, byrow = T))
    raster <- mask(raster, om.r)
    
    points <- points[vals > q, ]
  }
  
  r <- ZeroOneNorm(raster)
  r <- round(r, 2)

  thres <- seq(0, 1, len = 101)
  
  #Thresholding suitability layer
  r.thr <- r >= thres
  area.pred <- global(r.thr, mean, na.rm = T)$mean
  dArea <- area.pred[1:100] - area.pred[2:101]
  
  points.r <- as.data.frame(r.thr, xy = T)[, c("x", "y")]
  
  mp <- matrix(0, nrow = iterations, ncol = length(area.pred))
  mr <- matrix(0, nrow = iterations, ncol = length(area.pred))
  
  areas <- matrix(0, nrow = iterations, ncol = 4)
  colnames(areas) <- c("iteration", 
                       "Area.pres",
                       "Area.rand",
                       "PartialROC")
  
  for(i in 1:iterations){
    samp.pres <- sample(1:nrow(points), size = nrow(points) * p.points, replace = F)
    samp.rand <- sample(1:nrow(points.r), size = nrow(points) * p.points, replace = F)

    pres.thr <- terra::extract(r.thr, points[samp.pres, 1:2])[, -1]
    rand.thr <- terra::extract(r.thr, points.r[samp.rand, ])[, -1]
    
    
    #Calculating proportion of predicted points  
    means.pres <- colMeans(pres.thr)
    means.rand <- colMeans(rand.thr)
    
    #Calculating areas
    
    rects.pres <- dArea * means.pres[-101]
    rects.rand <- dArea * means.rand[-101]
    
    Area.pres <- sum(rects.pres)
    Area.rand <- sum(rects.rand)
    Area.ratio <- Area.pres / Area.rand
    
    areas[i, ] <- c(iteration = i, 
                    Area.pres = Area.pres,
                    Area.rand = Area.rand,
                    PartialROC = Area.ratio)
    mp[i, ] <- means.pres
    mr[i, ] <- means.rand
  }
  
  areas  <- as.data.frame(areas)

  rat <- round(mean(areas$PartialROC), 2)
  P <- with(areas, length(which(PartialROC < 1))/iterations)
  
  if(save.pdf){
  pdf(paste0(pdf.name, ".pdf"), width = 5, height = 5)
    plot(area.pred, colMeans(mp), xlab = "% Area predicted", 
         ylab = "1 - Omission error", col = "grey95", type = "l", 
         xlim = c(0, 1), ylim = c(0, 1), main = paste0("AUC ratio = ", rat, "\n P = ", P))
    lines(area.pred, colMeans(mr), col = "grey95", lty = 2, type = "l")
    for(j in 1:iterations){
      lines(area.pred, mp[j, ], col = "grey95", lwd = 0.25, type = "l")
      lines(area.pred, mr[j, ], col = "grey95", lty = 2, lwd = 0.25, type = "l")
    }
    lines(area.pred, colMeans(mp, na.rm = T), col = "red", lwd = 1.5, type = "s")
    lines(area.pred, colMeans(mr, na.rm = T), type = "s")
  dev.off()
}

  
  areas <- as.data.frame(areas)
  
  return(areas)
}
