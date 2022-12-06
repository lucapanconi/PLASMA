#Generalised Polarisation Calculations.
#' GP Calculate
#' 
#' Calculates the GP values of a point cloud given the photon counts.
#' 
#' @param points Point pattern including photon counts in columns 3 and 4.
#' @return Vector of GP values.
#' @export
gpcalculate <- function(points){
  #Calculate GP values.
  gpvalues <- (points[,3] - points[,4])/(points[,3] + points[,4])
  #Remove NA values.
  gpvalues[is.na(gpvalues)] <- 0
  #Return GP values.
  return(gpvalues)
}

#' GP Histogram
#' 
#' Creates a histogram of given GP values.
#' 
#' @param gpvalues A vector of GP values.
#' @return Histogram of GP values.
#' @export
gphistogram <- function(gpvalues, binwidth = 0.01, density = FALSE, saveplot = FALSE, location = ""){
  #Create dataframe.
  gpvalues <- as.data.frame(gpvalues)
  colnames(gpvalues) <- c("gp")
  p <- ggplot2::geom_histogram(binwidth = binwidth, color = "darkgrey", fill = "darkgrey")
  #Use density.
  if(density){
    #Create plot.
    ylab <- "Density"
    titl <- "Density of GP Values"
    print(
      ggplot2::ggplot(data = gpvalues, ggplot2::aes(x = gp, y = ..density..)) + p +
        ggplot2::geom_density(alpha = 0, fill = "#000000", size = 1) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = "GP Value", y = ylab) +
        ggplot2::ggtitle(titl) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                       axis.title.x = ggplot2::element_text(size = 16), axis.title.y = ggplot2::element_text(size = 16),
                       axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_text(size = 12),
                       plot.title = ggplot2::element_text(size = 20, hjust = 0.5))
    )
    #Save plot.
    if(saveplot){
      ggplot2::ggsave(paste(location, titl, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
    }
  } else {
    #Create plot.
    ylab <- "Count"
    titl <- "Histogram of GP Values"
    print(
      ggplot2::ggplot(data = gpvalues, ggplot2::aes(x = gp)) + p +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = "GP Value", y = ylab) +
        ggplot2::ggtitle(titl) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                       axis.title.x = ggplot2::element_text(size = 16), axis.title.y = ggplot2::element_text(size = 16),
                       axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_text(size = 12),
                       plot.title = ggplot2::element_text(size = 20, hjust = 0.5))
    )
    #Save plot.
    if(saveplot){
      ggplot2::ggsave(paste(location, titl, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
    }
  }
}

#' Ripley's H Function
#' 
#' Calculates Ripley’s H function using the spatstat Ripley’s K Function with isotropic edge correction and gives a line plot. 
#' 
#' @param points Point pattern.
#' @return Plot of Ripley's H Function.
#' @export
ripleyh <- function(points, window = NULL, createplot = TRUE, saveplot = FALSE, location = ""){
  #Determine save plot status.
  if(!createplot){
    saveplot <- FALSE
  }
  #Extract spatial coordinates only. If window is null, calculate window.
  if(is.null(window)){
    window <- c(min(points[,1]), max(points[,1]), min(points[,2]), max(points[,2]))
  }
  points <- spatstat.geom::ppp(points[,1],points[,2], c(window[1], window[2]), c(window[3], window[4]))
  #Get Ripley's K with isotropic edge correction.
  K <- as.matrix(spatstat.core::Kest(points, correction = c("isotropic"), nlarge = Inf))
  #Calculate Ripley's H.
  Hr <- as.data.frame(cbind(K[,1], sqrt(K[,3]/pi) - K[,1]))
  colnames(Hr) <- c("r", "H")
  #Plot.
  if(createplot){
    print(ggplot2::ggplot(data = Hr, ggplot2::aes(x = r)) +
            ggplot2::geom_line(ggplot2::aes(y = H), size = 1) +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
            ggplot2::labs(x = "r", y = "H(r)") +
            ggplot2::ggtitle("Ripley's H Function") +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5),
                           panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                           axis.line = ggplot2::element_line(colour = "black"), axis.title.x = ggplot2::element_text(size = 16),
                           axis.title.y = ggplot2::element_text(size = 16)))
    #Save plot.
    if(saveplot){
      ggplot2::ggsave(paste(location, "Ripley's H Function.pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
    }
  }
  #Return values.
  return(Hr)
}

#' Inform Radius
#' 
#' Determine an appropriate search radius based on the argmax of the Ripley’s H function.
#' 
#' @param points Point pattern.
#' @return Search radius.
#' @export
informradius <- function(points){
  #Get H values.
  Hr <- ripleyh(points, window = NULL, createplot = FALSE)
  #Find and return argmax.
  return(Hr[which.max(Hr[,2]),1])
}

#' Calculate Intersections
#' 
#' Determines point of intersection between two Gaussians.
#' 
#' @param m1 Mean of first Gaussian.
#' @param s1 Standard deviation of first Gaussian.
#' @param c1 Coefficient of first Gaussian.
#' @param m2 Mean of second Gaussian.
#' @param s2 Standard deviation of second Gaussian.
#' @param c2 Coefficient of second Gaussian.
#' @return Point of intersection.
#' @export
calcinter <- function(m1, m2, s1, s2, c1, c2){
  a <- 1/(2 * s2 ** 2) - 1/(2 * s1 ** 2)
  b <- m1/(s1 ** 2) - m2/(s2 ** 2)
  c <- log(c1/c2) + log(s2/s1) + (m2 ** 2)/(2 * s2 ** 2) - (m1 ** 2)/(2 * s1 ** 2) 
  return(pracma::roots(c(a, b, c)))
}

#' Inform Deviance
#' 
#' Try to determine an appropriate deviance threshold for JOSEPH based on the distribution of marked values.
#' 
#' @param values Marked values.
#' @return Deviance threshold.
#' @export
informdeviance <- function(values, k = 2, createplot = FALSE, saveplot = FALSE, location = "", outputall = FALSE, plotresolution = 100){
  #Determine save plot status.
  if(!createplot){
    saveplot <- FALSE
  }
  #Create Gaussian mixture model using Mclust and extract parameters.
  parameters <- mclust::Mclust(values, G = k)[[13]]
  #Extract mean.
  means <- as.numeric(do.call(c, lapply(parameters[[2]], function(i){
    return(as.numeric(i))
  })))
  #Check if means have been found.
  if(length(means) == 0){
    print("Parameter Estimation: No means found, standard deviation is likely 0.")
    return(NULL)
  } else {
    #Find other parameters.
    sds <- as.numeric(do.call(c, lapply(parameters[[3]][[4]], function(i){
      return(sqrt(as.numeric(i)))
    })))
    if(length(sds) < k){
      sds <- rep(sds[1], k)
    }
    coeff <- as.numeric(do.call(c, lapply(parameters[[1]], function(i){
      return(as.numeric(i))
    })))
    #Between each pair of consecutive curves, determine the points of intersection
    roots <- do.call(list, lapply(2:k, function(i){
      #Calculate all roots.
      allroots <- calcinter(means[i - 1], means[i], sds[i - 1], sds[i], coeff[i - 1], coeff[i])
      #Check cases.
      if(length(allroots) == 0){
        #No roots found, return NULL.
        print("Parameter Estimation: No roots found. This suggests data is not Normally distributed.")
        return(NULL)
      } else {
        if(is.complex(allroots) || sum(allroots >= means[i - 1] & allroots <= means[i]) == 0){
          #Roots are complex or outside of peaks, return mean position of peaks.
          root <-  mean(means[(i - 1):i])
          return(list(allroots, root, min(abs(root - means))))
        } else {
          #Calculate root between means.
          root <- allroots[allroots >= means[i - 1] & allroots <= means[i]]
          return(list(allroots, root, min(abs(root - means))))
        }
      }
    }))
    
    #Find minimum deviance and split point.
    deviances <- do.call(c, lapply(roots, function(i){
      if(is.null(i)){
        return(Inf)
      } else {
        return(i[[3]])
      }
    }))
    deviance <- min(deviances)
    split <- roots[[which.min(deviances)]][[2]]
    
    #Extract all intersections.
    intersection <- do.call(c, lapply(roots, function(i){
      if(is.null(i)){
        return(NULL)
      } else {
        return(i[[1]])
      }
    }))
    intersection <- intersection[!is.null(intersection)]
    
    #Create plot.
    if(createplot){
      #Create x data.
      x <- pracma::linspace(min(values), max(values), plotresolution)
      data <- data.frame()
      #Add Gaussians.
      for(i in 1:k){
        data <- rbind(data, cbind(x, paste("y", i, sep = ""), coeff[i] * dnorm(x, mean = means[i], sd = sds[i])))
      }
      #Update column names.
      colnames(data) <- c("gp", "var", "dens")
      #Return to numeric.
      data[,1] <- as.numeric(data[,1])
      data[,3] <- as.numeric(data[,3])
      #Create plot.
      print(
        p <- ggplot2::ggplot(data, ggplot2::aes(x = gp, y = dens)) +
          ggplot2::geom_line(ggplot2::aes(color = var), size = 1) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::labs(x = "GP Value", y = "Density") +
          ggplot2::ggtitle("Estimated Distributions of GP Values") +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5),
                         panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(colour = "black"), axis.title.x = ggplot2::element_text(size = 16),
                         axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_text(size = 12),
                         axis.title.y = ggplot2::element_text(size = 16), legend.position = "none")
      )
      #Save plot.
      if(saveplot){
        ggplot2::ggsave(paste(location, "Estimated Distributions of GP Values.pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
      }
    }
    
    #Output data.
    if(outputall){
      return(list(mean = means, sd = sds, coefficient = coeff, intersection = intersection, split = split, deviance = deviance))
    } else {
      return(deviance)
    }
  }
}