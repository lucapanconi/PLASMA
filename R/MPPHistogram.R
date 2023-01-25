#Marked Point Pattern Histogram.
#Bin Datum function used in MPP Histogram
bindatum <- function(i, y_wid, GP_x, GP_y, pts){
  return(c(1 + (i - 1) %/% y_wid, 1 + (i - 1) %% y_wid,
           sum(pts[GP_x == 1 + (i - 1) %/% y_wid & GP_y == 1 + (i - 1) %% y_wid, 3]),
           sum(pts[GP_x == 1 + (i - 1) %/% y_wid & GP_y == 1 + (i - 1) %% y_wid, 4])))
}

#Compile Matrix function used in MPP Histogram.
compilematrix <- function(i, rep, col){
  return(rep[rep[,1] == i, col])
}

#MPP Histogram.
mom <- function(points, binsize = 50, location = "", plotname = "GP Histogram", useggplot = FALSE, outputdata = FALSE){
  #If GP values is not given, calculate it.
  if(ncol(points) == 4){
    points <- cbind(points, gpcalculate(points))
  }
  
  #Calculate minimum and maximum x and y values.
  min_x <- min(points[,1])
  min_y <- min(points[,2])
  max_x <- max(points[,1])
  max_y <- max(points[,2])
  
  #Update column names.
  colnames(points) <- c("x", "y", "green_photons", "red_photons", "GP")
  
  #Shift coordinates to the origin.
  points[,1] <- points[,1] - min_x
  points[,2] <- points[,2] - min_y
  
  #Calculate size of 2D histogram.
  x_width <- ceiling((max_x - min_x) / binsize)
  y_width <- ceiling((max_y - min_y) / binsize)
  
  #Divide each x and y coordinate by the bin size and round up.
  combGP_x <- ceiling(points[,1] / binsize)
  combGP_y <- ceiling(points[,2] / binsize)
  
  #Put zeroes into first entry.
  combGP_x[combGP_x == 0] <- 1
  combGP_y[combGP_y == 0] <- 1
  
  #Define GP histogram bins.
  GP_histbins <- seq(-1,1,0.05)
  
  #Create images.
  if(useggplot){
    #Create data frame for individual bins.
    bindata <- as.data.frame(do.call(rbind, lapply(1:(x_width * y_width), bindatum, y_wid = y_width, GP_x = combGP_x, GP_y = combGP_y, pts = points)))
    
    #Add GP.
    gpvalues <- (bindata[,3] - bindata[,4])/(bindata[,3] + bindata[,4])
    bindata <- cbind(bindata, gpvalues)
    
    #Update 0 photon regions.
    bindata[bindata[,3] == 0, 3] <- NA
    bindata[bindata[,4] == 0, 4] <- NA
    
    #Update column names.
    colnames(bindata) <- c("x", "y", "channel1", "channel2", "GP")
    
    #Create GP image matrix.
    GP_image_cols <- lapply(1:x_width, compilematrix, rep = bindata, col = 5)
    GP_image <- t(do.call(cbind, GP_image_cols))
    
    #Output data or plot.
    if(outputdata){
      #Output bin data.
      return(list(bindata, GP_image))
    } else {
      #Create with ggplot.
      pdf(paste(location, plotname, ".pdf", sep = ""), onefile = TRUE, width = 8)
      hist(points$green_photons, xlab = "Photon Count", ylab = "Frequency", main = "Histogram of Photon Count in First Channel")
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
      hist(points$red_photons, xlab = "Photon Count", ylab = "Frequency", main = "Histogram of Photon Count in Second Channel")
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
  } else {
    #Create photon matrices.
    GP_gph_im <- matrix(0, x_width, y_width)
    GP_rph_im <- matrix(0, x_width, y_width)
    for(i in 1:nrow(points)){
      GP_gph_im[combGP_x[i], combGP_y[i]] <- GP_gph_im[combGP_x[i], combGP_y[i]] + points[i,3]
      GP_rph_im[combGP_x[i], combGP_y[i]] <- GP_rph_im[combGP_x[i], combGP_y[i]] + points[i,4]
    }

    #Replace zeroes.
    GP_gph_im[GP_gph_im == 0] <- NA
    GP_rph_im[GP_rph_im == 0] <- NA

    #Create images..
    GP_image <- (GP_gph_im - GP_rph_im)/(GP_gph_im + GP_rph_im)
    #Create plot.
    pdf(paste(location, plotname, ".pdf", sep = ""))
    #Create histograms.
    hist(points$green_photons)
    fields::image.plot(GP_gph_im)
    hist(points$red_photons)
    fields::image.plot(GP_rph_im)
    hist(GP_image, GP_histbins)
    gpcol <- colorRampPalette(c("cyan","magenta"))
    GPcol <- gpcol(10)[cut(c(-1,1), breaks = 10)]
    fields::image.plot(GP_image,col = gpcol(10))
    dev.off()
  }
  #Save important data.
  save.image(paste(location, plotname, ".Rdata", sep = ""))
  MASS::write.matrix(GP_image, paste(location, plotname, ".txt", sep = ""))
}