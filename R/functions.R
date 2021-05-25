# all the functions

#### Scales
# Determining scale from files

# Function to calculate scale from 2 points (in pixels) and known distance:
# modified version to get the coordinates directly from the scale pos file. Input, the file name and the known distance.

CalcPixelScale <- function(sc.coords, distance) {
  # Pixel distance between the two points
  pixels <- sqrt((sc.coords$X.1 - sc.coords$X)^2 + (sc.coords$Y.1 - sc.coords$Y)^2)
  # Return scale: real world distance divided by pixel distance
  distance / pixels
}


#### Load samples
# Loads both the scale and pos file, fixes orientation of the coordinates, creates a trajectory file from the coords, creates the scale from the scale file, it scales the trajectory object.

# ```{r Sample from files}

# loadSample <- function (sampleFile, scaleFile, knownScale) {
#   # get coords file
#   coords <- read.table(sampleFile, skip = 2, header = T)
#   # get the scale file
#   sc.coords <- read.table(scaleFile, sep = "", skip = 2, header = T)
#   # invert the direction of Y column values = y-axis so the trajectory is not upside-down
#   coords$Y.3 <- max(coords$Y.3) - coords$Y.3
#   
#   # create trajectory from coordinates
#   # columns as included in the pos files direct from MathLab
#   trj <- TrajFromCoords(coords, xCol = 8, yCol = 9, fps = 250, spatialUnits = "pixels")
#   # run function, relative to the desired unit to use
#   # example, known d = 17mm, to meters then: 0.017, to cm then: 1.7
#   scale <- CalcPixelScale(sc.coords, knownScale)
#   TrajScale(trj, scale, "m")
# }


# Plots acceleration and speed together - not used in this case, see next.

# plotAcSp <- function(derivs) {
#   # Plot acceleration and speed
#   plot(derivs$acceleration ~ derivs$accelerationTimes, type = 'l', col = 'red', 
#        yaxt = 'n',
#        xlab = 'Time (s)',
#        ylab = expression(paste('Acceleration (', m/s^2, ')')))
#   axis(side = 2, col = "red")
#   lines(derivs$speed ~ derivs$speedTimes, col = 'blue')
#   axis(side = 4, col = "blue")
#   mtext('Speed (m/s)', side = 4, line = 3)
#   abline(h = 0, col = 'lightGrey')
# }




## Load from List
# Note that fill = T has been added to take in the files from chronos cam
# Needs to be checked for the files from optronix

loadSample_FromList <- function (listObject, number) {
  
  # get coords file
  coords <- read.table(listObject[[1]][number], skip = 2, header = T, fill = T)
  coords$filename <- listObject[[1]][number]
  
  # get the scale file
  sc.coords <- read.table(listObject[[2]][number], sep = "", skip = 2, header = T)
  sc.coords$filename <- listObject[[2]][number]
  
  # invert the direction of Y column values = y-axis so the trajectory is not upside-down
  coords$Y.3 <- max(coords$Y.3) - coords$Y.3
  
  # create trajectory from coordinates
  # columns as included in the pos files direct from MathLab
  trj <- TrajFromCoords(coords[, c(8, 9)], xCol = 1, yCol = 2, fps = 250, spatialUnits = "pixels")
  # run function, relative to the desired unit to use
  # example, known d = 17mm, to meters then: 0.017, to cm then: 1.7
  scale <- CalcPixelScale(sc.coords, listObject[[3]]$scaleValues[number])
  TrajScale(trj, scale, "m")
}

loadSample_FromList2 <- function (listObject1, listObject2, listObject3) {
  
  # get coords file
  coords <- read.table(listObject1, skip = 2, header = T, fill = T)
  # coords$filename <- listObject1
  
  # get the scale file
  sc.coords <- read.table(listObject2, sep = "", skip = 2, header = T)
  # sc.coords$filename <- listObject2
  
  # invert the direction of Y column values = y-axis so the trajectory is not upside-down
  coords$Y.3 <- max(coords$Y.3) - coords$Y.3
  
  # create trajectory from coordinates
  # columns as included in the pos files direct from MathLab
  trj <- TrajFromCoords(coords[, c(8, 9)], xCol = 1, yCol = 2, fps = 250, spatialUnits = "pixels")
  # run function, relative to the desired unit to use
  # example, known d = 17mm, to meters then: 0.017, to cm then: 1.7
  scale <- CalcPixelScale(sc.coords, listObject3)
  TrajScale(trj, scale, "m")
}



# sequential labels with prefix from 0 to 99
mylabelseq <- function(prefix, from, to, by) {
  x <- seq(from, to, by)
  paste0(prefix, ifelse(str_detect(x, "[0-9]{2}"), x, paste0("0",x)))
}


RightOrNot <- function(trajectoryO, title) {
  # par(mfrow = c(1, 1))
  
  smoothedO <- TrajSmoothSG(trajectoryO, p = 3, n = 21)
  smootrjO <- TrajDerivatives(smoothedO)
  
  # the prep
  startIndex_O <- which(smootrjO$speed > 1e-8)[1] - 1 
  endIndex_O <- startIndex_O + (150/4)
  local_maxima_text_O <- JLocalMaxima(smootrjO$speed, 20/4, startIndex_O, endIndex_O)
  
  # the check
  plot(smootrjO$speed, type = "l")
  title(title)
  abline(v = local_maxima_text_O, col = "red")
  abline(v = c(startIndex_O, endIndex_O), col = "blue")
  
}

## Right from Me


# smooth, derivatives, plot speed
RightFromMe <- function(trajectoryO, firstMoveFrame, contact, title) {
  # par(mfrow = c(1, 1))
  
  smoothedO <- TrajSmoothSG(trajectoryO, p = 3, n = 21)
  smootrjO <- TrajDerivatives(smoothedO)
  
  # the prep
  startIndex_O <- firstMoveFrame 
  endIndex_O <- startIndex_O + (150/4)
  local_maxima_text_O <- JLocalMaxima(smootrjO$speed, 20/4, startIndex_O, endIndex_O)
  
  # the check
  plot(smootrjO$speed, type = "l")
  title(title)
  abline(v = local_maxima_text_O, col = "red")
  abline(v = c(startIndex_O, endIndex_O), col = "blue")
  abline(v = startIndex_O + contact, col = "green")
  
}



## Beautifully pure

def_smoo_set <- c(3, 21)
alt_smoo_set <- c(4, 23)

def_start_threshold <- 1e-8
alt_start_threshold <- 1e-3



DerivsOnly <- function(trajectoryO, smoothSettings) {
  
  smoothedO <- TrajSmoothSG(trajectoryO, p = smoothSettings[1], n = smoothSettings[2])
  smootrjO <- TrajDerivatives(smoothedO) 
  smootrjO
}


IndexLimitsAndPlot <- function(derivs_object, startThreshold = 1e-8, limit_range = 150, window_size = 20) {
  # the prep
  startIndex_O <- which(derivs_object$speed > startThreshold)[1] - 1 
  endIndex_O <- startIndex_O + (limit_range/4)
  local_maxima_text_O <- as.vector(JLocalMaxima(derivs_object$speed, window_size/4, startIndex_O, endIndex_O))
  
  maximaIndexes <- list(startIndex_O, endIndex_O, local_maxima_text_O)
  names(maximaIndexes) <- c("start_index", "end_index", "maxim(um/a)")
  
  return(maximaIndexes)
  
}


IndexLimitsManual <- function(derivs_object, startThreshold = 1e-8, thresholdFrame = 1, limit_range = 155, window_size = 35) {
  # the prep
  startIndex_O <- which(derivs_object$speed > startThreshold)[thresholdFrame] - 1 
  endIndex_O <- startIndex_O + (limit_range/4)
  local_maxima_text_O <- as.vector(JLocalMaxima(derivs_object$speed, window_size/4, startIndex_O, endIndex_O))
  
  maximaIndexes <- list(startIndex_O, endIndex_O, local_maxima_text_O)
  names(maximaIndexes) <- c("start_index", "end_index", "maxim(um/a)")
  
  return(maximaIndexes)
  
}


DerivsAndPlot <- function(trajectoryO, smoothSettings, startThreshold = 1e-8, title, plot = FALSE, trim = FALSE) {
  
  smoothedO <- TrajSmoothSG(trajectoryO, p = smoothSettings[1], n = smoothSettings[2])
  smootrjO <- TrajDerivatives(smoothedO) 
  
  # the prep
  startIndex_O <- which(smootrjO$speed > startThreshold)[1] - 1 
  endIndex_O <- startIndex_O + (150/4)
  local_maxima_text_O <- JLocalMaxima(smootrjO$speed, 20/4, startIndex_O, endIndex_O)
  
  if (plot == FALSE) {
    
    maximaIndexes <- c(startIndex_O, endIndex_O, local_maxima_text_O)
    names(maximaIndexes) <- c("start_index", "end_index", "maximum")
    
    outcome1 <- list(smootrjO, maximaIndexes)
    
    return(outcome1)
    
  }
  else {
    if (trim == FALSE) {
      # the plot to check
      plot(smootrjO$speed, type = "l")
      title(title)
      abline(v = local_maxima_text_O, col = "red")
      abline(v = c(startIndex_O, endIndex_O), col = "blue")
      
    }
    
    else {
      # trim = TRUE
      plot(smootrjO$speed[startIndex_O:local_maxima_text_O], type = "l", 
           xlim = c(0, 125))
      title(title)
      abline(v = local_maxima_text_O, col = "red")
      abline(v = c(startIndex_O, endIndex_O), col = "blue")
      
    }
    
  }
}


# allows manual input of start index, end index (length after start)
DerivsAndPlot_manual <- function(trajectoryO, smoothSettings, start_index, end_index, title, plot = FALSE) {
  
  smoothedO <- TrajSmoothSG(trajectoryO, p = smoothSettings[1], n = smoothSettings[2])
  smootrjO <- TrajDerivatives(smoothedO) 
  
  # the prep
  startIndex_O <- start_index 
  endIndex_O <- startIndex_O + (end_index/4)
  local_maxima_text_O <- JLocalMaxima(smootrjO$speed, 20/4, startIndex_O, endIndex_O)
  
  if (plot == FALSE) {
    
    maximaIndexes <- c(startIndex_O, endIndex_O, local_maxima_text_O)
    names(maximaIndexes) <- c("start_index", "end_index", "maximum")
    
    outcome1 <- list(smootrjO, maximaIndexes)
    
    return(outcome1)
    
  }
  else {
    
    # the plot to check
    plot(smootrjO$speed, type = "l")
    title(title)
    abline(v = local_maxima_text_O, col = "red")
    abline(v = c(startIndex_O, endIndex_O), col = "blue")
    
  }
}



#### Speed plot
# Plot for the speed change during the trajectory
# Need the derivs (derivatives) files to run
# only speed 
# function to plot only speed
plotSpeed <- function(derivs) {
  plot(derivs$speed ~ derivs$speedTimes, type = 'l', col = 'black', 
       yaxt = 'n',
       xaxt = "n", 
       xlab = 'Time (s)',
       ylab = expression(paste('Speed (', m/s, ')')))
  axis(side = 2, col = "black")
  axis(side = 1, at=c(seq(from = min(derivs$speedTimes),to = max(derivs$speedTimes), by = 0.1)))
}



#### Local Maxima

JLocalMaxima <- function(v, window = 1, startIndex = 1, endIndex = length(v))
{
  getWindow <- function(i) {
    # Don't try to look past the ends of the data
    si <- max(1, i - window)
    ei <- min(length(v), i + window)
    v[si : ei]
  }
  
  maxima <- numeric(length(v) / 2)
  nm <- 0
  for (i in startIndex:endIndex) {
    # Is this point a maximum?
    if (v[i] == max(getWindow(i))) {
      nm <- nm + 1
      maxima[nm] <- i
    }
  }
  
  utils::head(maxima, nm)
}



## the RightOrNot

# smooth, derivatives, plot speed

SpeedAndTime2 <- function(derivatives_list, multiMax = 1, ...) {
  
  indexes <- IndexLimitsManual(derivs_object = derivatives_list, ...)
  
  max_speed <- (derivatives_list$speed[indexes$`maxim(um/a)`[multiMax]])*100
  
  elapsed_time <- (derivatives_list$speedTimes[indexes$`maxim(um/a)`[multiMax]] - derivatives_list$speedTimes[indexes$start_index])*1000
  
  tibble(max_speed, elapsed_time)
  
}

R2xtCorrCoeff <- function(lvar, cvar) {
  rxc <- cor(lvar, cos(cvar)) ; rxs <- cor(lvar, sin(cvar))
  rcs <- cor(cos(cvar), sin(cvar))
  R2xtVal <- ((rxc*rxc)+(rxs*rxs)-(2*rxc*rxs*rcs))/(1-rcs*rcs)
  return(R2xtVal)
}

R2xtIndTestRand <- function(lvar, cvar, NR) {
  R2xtObs <- R2xtCorrCoeff(lvar, cvar); nxtrm <- 1
  for (r in 1:NR) {
    lvarRand <- sample(lvar)
    R2xtRand <- R2xtCorrCoeff(lvarRand, cvar) 
    if (R2xtRand >= R2xtObs) { nxtrm <- nxtrm+1 } }
  pval <- nxtrm/(NR+1)
  names(R2xtObs) <- "coefficient"
  names(pval) <- "p-value"
  return(c(R2xtObs, pval))
}
