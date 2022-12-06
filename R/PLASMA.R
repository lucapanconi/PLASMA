#All PLASMA Functions.
#' PLASMA
#' 
#' Performs P-Check, JOSEPH and MPP Histogram on a point pattern or list of point patterns.
#' 
#' @param points Point pattern or list of point patterns.
#' @param location Output location for results.
#' @return Results of P-Check, JOSEPH and MPP Histogram, including: histogram of GP values with density plot, density plots of fitted Gaussians with all Gaussian parameters, Ripley’s H Function, all suggested parameters, all P-Check p-values, all JOSEPH results (points data), a plot of all JOSEPH clusters, summary statistics of clusters found by JOSEPH over all clouds (points per cluster, area of convex/concave hull of cluster, density, mean mark per cluster, variance of mark values per cluster, global mean mark, difference between global mean and cluster mean), all MPP Histogram results.

#' @export
plasma <- function(points, location = "", radius = NULL, deviance = NULL, split = NULL, minpoints = 2, names = NULL){
  #Convert single cloud to list.
  if(!is.null(ncol(points))){
    points <- list(points)
  }
  
  #Create names if not given.
  if(is.null(names)){
    names <- 1:length(points)
  }
  
  #Check if GP values need to be calculated.
  points <- lapply(points, function(cloud){
    if(ncol(cloud) == 4){
      #Calculate GP values from photon channels and add to data frame.
      originalnames <- colnames(cloud)
      cloud <- cbind(cloud, gpcalculate(cloud))
      colnames(cloud) <- c(originalnames, "GP")
      return(cloud)
    } else if(ncol(cloud) != 5){
      stop("Input data frame with four columns: x coordinates, y coordinates, photon count first channel and photon count second channel.")
    } else {
      return(cloud)
    }
  })
  
  #Compile all GP values.
  gpvalues <- do.call(c, lapply(points, function(cloud){
    return(gpcalculate(cloud))
  }))
  
  #Plot calculated GP values as a histogram/density plot.
  gphistogram(gpvalues = gpvalues, binwidth = 0.01, density = TRUE, saveplot = TRUE, location = location)
  
  #Split the GP distribution via mixture modelling and suggest a sensible deviance for JOSEPH.
  if(is.null(deviance)){
    parameters <- informdeviance(gpvalues, createplot = TRUE, saveplot = TRUE, location = location, outputall = TRUE)
    deviance <- as.numeric(parameters[[6]])
    saveRDS(parameters, paste(location, "Fitted Gaussian Parameters.Rdata", sep = ""))
  }
  
  #Analyse each cloud with each algorithm.
  bulkresults <- lapply(points, bulkanalyse, wholeradius = radius, wholedeviance = deviance, split = split, minpoints = minpoints)
  saveRDS(bulkresults, paste(location, "All Results.RData", sep = ""))
  
  #Get JOSEPH cluster statistics.
  josephresults <- lapply(bulkresults, extractresults, resnum = 2)
  josephclusterstatistics <- do.call(rbind, lapply(1:length(josephresults), function(i){
    res <- josephclustersummary(josephresults[[i]])
    if(!is.null(res)){
      return(cbind(names[i], res))
    }}))
  if(ncol(josephclusterstatistics) > 1){
    colnames(josephclusterstatistics) <- c("name", "number of points", "area", "density", "average mark", "mark standard deviation", "global average", "difference")
    josephclusterstatistics <- josephclusterstatistics[!is.null(josephclusterstatistics[,2]),]
    write.csv(josephclusterstatistics, paste(location, "JOSEPH Cluster Summary Data.csv", sep = ""), row.names = FALSE)
  }
  
  #Plot JOSEPH results.
  pdf(paste(location, "JOSEPH Cluster Results.pdf", sep = ""), onefile = TRUE)
  for(jclust in josephresults){
    if(!is.character(jclust)){
      print(RSMLM::plotClusterScatter(jclust[,1:2], jclust[,4]))
    }
  }
  dev.off()
  
  #Plot JOSEPH Hulls.
  pdf(paste(location, "JOSEPH Cluster Hulls.pdf", sep = ""), onefile = TRUE)
  for(jclust in josephresults){
    if(!is.character(jclust)){
      minx <- min(jclust[,1])
      maxx <- max(jclust[,1])
      miny <- min(jclust[,2])
      maxy <- max(jclust[,2])
      josephvers <- data.frame()
      if(length(unique(jclust[,4])) > 1){
        for(i in unique(jclust[,4])){
          if(i > 0){
            joe <- jclust[jclust[,4] == i,1:2]
            josephvers <- rbind(josephvers, cbind(joe[chull(joe),], i))
          }
        }
        colnames(josephvers) <- c("x", "y", "v")
        print(ggplot2::ggplot(josephvers, ggplot2::aes(x = x, y = y, fill = factor(v))) +
                ggplot2::geom_polygon() +
                ggplot2::scale_x_continuous(limits = c(minx, maxx), expand = c(0.01, 0.01)) +
                ggplot2::scale_y_continuous(limits = c(miny, maxy), expand = c(0.01, 0.01)) +
                ggplot2::scale_fill_manual(values = rep("purple", max(josephvers[,3]))) +
                ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1),
                               axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                               axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                               axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                               legend.position = "none"))
      }
    }
  }
  dev.off()
  
  #Define GP histogram bins.
  histbins <- seq(-1,1,0.05)
  
  #Create PDF with all MPP Histogram results and isolate GP images.
  mpphistogramresults <- lapply(bulkresults, extractresults, resnum = 3)
  lapply(1:length(mpphistogramresults), plothistogramresults, bindata = mpphistogramresults, GP_histbins = histbins, location = location, points = points)
  
  #Save GP images.
  saveRDS(lapply(mpphistogramresults, function(bindata){
    if(!is.character(bindata)){
      return(bindata[[2]])
    }
  }), paste(location, "GP Images.RData", sep = ""))
}

#' Bulk Analyse
#' 
#' Used in PLASMA.
#' 
#' @param cloud Point pattern or list of point patterns.
#' @return One set of results from P-Check, JOSEPH and MPP Histogram.
#' @export
bulkanalyse <- function(cloud, wholeradius, wholedeviance, split, minpoints){
  #Convert continuous point cloud into discretised point cloud.
  discretecloud <- discretise(cloud[,c(1,2,5)], splits = split, names = c("Ordered", "Disordered"))
  
  #Suggest search radius.
  if(is.null(wholeradius)){
    wholeradius <- informradius(discretecloud)
  }
  
  #Perform P-Check.
  presults <- pcheck(discretecloud, radius = wholeradius, outputall = TRUE)
  
  #If P-Check gives significant p-value, proceed with JOSEPH, otherwise do not.
  if(as.numeric(presults[[1]][2]) <= 0.05){
    #Perform JOSEPH on continuous point pattern.
    josephresults <- joseph(cloud[,c(1,2,5)], deviance = wholedeviance, radius = wholeradius)
    josephresults <- filterjoseph(josephresults = josephresults, mode = c("dispersion", "minpoints"), minpoints = minpoints)
    
    #Perform MPP Histogram.
    mpphistogramresults <- mpphistogram(cloud, location = location, plotname = paste("GP Histogram ", cloudnumber, sep = ""), useggplot = TRUE, outputdata = TRUE)
  } else {
    #Update data for JOSEPH and MPP Histogram results.
    josephresults <- paste("No evidence of separation found (p = ", presults[[1]][2], ")", sep = "")
    mpphistogramresults <- paste("No evidence of separation found (p = ", presults[[1]][2], ")", sep = "")
  }
  
  #Return results.
  return(list("P-Check Results" = presults, "Joseph Results" = josephresults, "MPP Histogram Results" = mpphistogramresults))
}

#' Extract Results.
#' 
#' Used in PLASMA.
#' 
#' @param res Results from Bulk Analyse.
#' @return One of the results from P-Check, JOSEPH or MPP Histogram.
#' @export
extractresults <- function(res, resnum){
  return(res[[resnum]])
}

#' Plot Histogram Results
#' 
#' Used in PLASMA.
#' 
#' @param i Index.
#' @return Plots of results from MPP Histogram.
#' @export
plothistogramresults <- function(i, bindata, GP_histbins, location, points){
  if(!is.character(bindata[[i]])){
    bindata <- bindata[[i]][[1]]
    cloud <- as.data.frame(points[[i]])
    colnames(cloud) <- c("x", "y", "channel1", "channel2", "GP")
    plotname <- paste("Plot ", i, sep = "")
    pdf(paste(location, "GP Histograms ", plotname, ".pdf", sep = ""), onefile = TRUE, width = 8)
    hist(cloud$channel1, xlab = "Photon Count", ylab = "Frequency", main = "Histogram of Photon Count in First Channel")
    print(
      ggplot2::ggplot(bindata, ggplot2::aes(x = x, y = y)) + ggplot2::geom_raster(ggplot2::aes(fill = channel1)) +
        ggplot2::scale_fill_viridis_c(name = "Photons", na.value = "white") +
        ggplot2::scale_x_continuous(expand = c(0,0)) + 
        ggplot2::scale_y_continuous(expand = c(0,0)) +
        ggplot2::xlab("x") + ggplot2::ylab("y") +
        ggplot2::ggtitle(paste("Number of Photons in First Channel across ROI (", plotname, ")", sep = "")) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                       axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14),
                       axis.text.x = ggplot2::element_text(size = 10), axis.text.y = ggplot2::element_text(size = 10),
                       legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.75, "cm"),
                       legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16))
    )
    hist(cloud$channel2, xlab = "Photon Count", ylab = "Frequency", main = "Histogram of Photon Count in Second Channel")
    print(
      ggplot2::ggplot(bindata, ggplot2::aes(x = x, y = y)) + ggplot2::geom_raster(ggplot2::aes(fill = channel2)) +
        ggplot2::scale_fill_viridis_c(name = "Photons", na.value = "white") +
        ggplot2::scale_x_continuous(expand = c(0,0)) + 
        ggplot2::scale_y_continuous(expand = c(0,0)) +
        ggplot2::xlab("x") + ggplot2::ylab("y") +
        ggplot2::ggtitle(paste("Number of Photons in Second Channel across ROI (", plotname, ")", sep = "")) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                       axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14),
                       axis.text.x = ggplot2::element_text(size = 10), axis.text.y = ggplot2::element_text(size = 10),
                       legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.75, "cm"),
                       legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16))
    )
    hist(bindata[,5], GP_histbins, xlab = "GP", ylab = "Frequency", main = "Histogram of GP Values")
    print(
      ggplot2::ggplot(bindata, ggplot2::aes(x = x, y = y)) + ggplot2::geom_raster(ggplot2::aes(fill = GP)) +
        ggplot2::scale_fill_viridis_c(name = "GP", na.value = "white") +
        ggplot2::scale_x_continuous(expand = c(0,0)) + 
        ggplot2::scale_y_continuous(expand = c(0,0)) +
        ggplot2::xlab("x") + ggplot2::ylab("y") +
        ggplot2::ggtitle(paste("GP Values across ROI (", plotname, ")", sep = "")) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                       axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14),
                       axis.text.x = ggplot2::element_text(size = 10), axis.text.y = ggplot2::element_text(size = 10),
                       legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.75, "cm"),
                       legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16))
    )
    dev.off()
  }
}