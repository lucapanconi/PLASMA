#All PLASMA Functions
require(mclust)
require(dplyr)

#Adjusted Rand Index Analysis.

#' Discretise
#'
#' Converts a continuous marked point pattern into a discrete marked point pattern.
#'
#' @param points A marked point pattern. Matrix or data frame with first three columns corresponding to x, y and mark.
#' @param splits Split values which define intervals to categorise points. If left as null, these intervals will be automatically estimated from Gaussian fitting.
#' @param names A vector of characters, where each entry corresponds to the name of a category.
#' @param k Number of Gaussians to fit in case splits is left null. Should be equal to number of suspected categories.
#' @param createplot Boolean dictating whether to create a plot of the new discretised marked point pattern.
#' @return A discretised marked point pattern.
#' @export
discretise <- function(points, splits = NULL, names = NULL, k = 2, createplot = TRUE){
  #Check if points has too many columns.
  points <- points[,1:3]
  #Split data set based on values in third column. If splits are not specified, guess using inform deviance.
  if(is.null(splits)){
    splits <- informdeviance(points[,3], k = k, createplot = FALSE, outputall = TRUE)[[5]]
    if(is.null(splits)){
      #Standard deviation is probably 0.
      uniq <- unique(points[,3])
      splits <- c()
      for(i in 2:length(uniq)){
        splits <- append(splits, mean(c(uniq[i-1], uniq[i])))
      }
    }
  }
  if(is.null(names)){
    names <- paste("Category ", 1:(length(splits) + 1), sep = "")
  }
  #Split data set.
  if(length(splits) == 1){
    thresh <- points[,3] <= as.numeric(splits[1])
    points[thresh,3] <- names[1]
    points[!thresh,3] <- names[2]
  } else {
    vals <- rep(0, nrow(points))
    vals[points[,3] <= as.numeric(splits[1])] <- names[1]
    for(i in 1:(length(splits) - 1)){
      vals[points[,3] > as.numeric(splits[i]) & points[,3] <= as.numeric(splits[i + 1])] <- names[i + 1]
    }
    vals[points[,3] > as.numeric(splits[length(splits)])] <- names[length(names)]
    points[,3] <- vals
  }
  #Update column names.
  colnames(points) <- c("x", "y", "category")
  #Create plot.
  if(createplot){
    print(
      ggplot2::ggplot(data = points, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(ggplot2::aes(color = factor(category))) +
        ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
        ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
        ggplot2::labs(color = "Category") +
        ggplot2::scale_colour_brewer(palette = "Set2") +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                       axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size = 20), legend.text = ggplot2::element_text(size = 16),
                       legend.position = "none")
    )
  }
  #Return new points.
  return(points)
}

#' Inform Cluster Persistence
#'
#' Estimates the persistence threshold to use for ToMATo clustering.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y.
#' @param radius Suspected cluster radius. If left null, this will be estimated from the Ripley's function.
#' @return A persistence threshold.
#' @export
informclusterpersistence <- function(points, radius = NULL){
  #Calculate number of neighbours within radius.
  if(is.null(radius)){
    radius <- informradius(points)
  }
  #Calculate number of neighbours within radius for each point.
  numberneighbours <- as.numeric(colSums(Rfast::Dist(points[,1:2]) <= radius))
  if(max(numberneighbours) < 2){
    print("Parameter Estimation: No neighbours found. This suggests that the given search radius is too small.")
    return(NULL)
  }
  #Calculate persistence.
  persistence <- informdeviance(numberneighbours[numberneighbours > 0])
  if(is.null(persistence)){
    return(NULL)
  } else {
    return(ceiling(persistence))
  }
}

#' Get ARI (Proxy for ToMATo Clustering)
#'
#' Perform ToMATo clustering.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y.
#' @param search_radius Cluster radius.
#' @param persistence Persistence threshold.
#' @return Point cloud with cluster index.
#' @export
getari <- function(points, search_radius, persistence){
  #Define useful properties.
  number_of_points <- nrow(points)
  distance_matrix <- as.matrix(Rfast::Dist(points[,1:2]))
  distance_radius <- distance_matrix <= search_radius
  #Get densities.
  densities <- as.numeric(colSums(distance_radius))
  #Get filtration.
  filtration <- 1:number_of_points
  point_ids <- Rfast::sort_cor_vectors(filtration, densities, descending = TRUE)
  filtration_ids <- Rfast::sort_cor_vectors(filtration, point_ids)
  #Perform persistent homology.
  r <- rep(0, number_of_points)
  for(i in filtration){
    #Get original point id.
    point_id <- point_ids[i]
    #Extract filtrations of neighbours of i.
    neighbour_filtrations <- filtration_ids[distance_radius[, point_id]]
    #Extract neighbours with lower filtration values than i.
    N <- neighbour_filtrations[neighbour_filtrations < i]
    #Check length of N.
    if(length(N) == 0){
      #If N is empty, then i is a peak. Create a new root. Update r and e.
      r[point_id] <- point_id
    }else{
      #If N is non-empty, then there must be at least one neighbour of higher density. #Cycle through neighbours and check persistence holds.
      potentials <- N[densities[point_ids[N]] <= densities[point_id] + persistence]
      #Check length of potentials.
      if(length(potentials) == 0){
        #Node cannot connect to any neighbours (persistence too low). Create new set for node.
        r[point_id] <- point_id
      } else {
        #Node can connect to at least one neighbour and bridges gap between other neighbours.
        potentials <- point_ids[Rfast::Sort(potentials)]
        j <- potentials[1]
        r[point_id] <- j
        if(length(potentials) > 1){
          #Merge the sets of the potentials.
          for(k in potentials[-1]){
            r_k <- r == k
            r[r_k] <- j
          }
        }
      }
    }
  }
  #Determine cluster indices.
  found <- 0
  #Filter out clusters with too few points.
  cluster_index <- rep(0, number_of_points)
  for(i in unique(r)){
    #Get points corresponding to that root.
    cluster_points <- filtration[r == i]
    #Filter out clusters below persistence.
    if(densities[i] <= persistence){
      #Set to outliers.
      cluster_index[cluster_points] <- 0
    } else {
      #Initialise new cluster.
      found <- found + 1
      cluster_index[cluster_points] <- found
    }
  }
  #Return new point cloud with cluster index.
  return(cbind(points, cluster_index))
}

#Generalised Polarisation Calculations.

#' GP Calculate
#'
#' Calculate the GP values from a marked point pattern.
#'
#' @param points A marked point pattern. Matrix or data frame with third and fourth columns corresponding to green and red photon counts.
#' @return Vector of GP values of each localisation.
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
#' Create histogram of all GP values.
#'
#' @param gpvalues A vector of GP values.
#' @param binwidth Width of bin used in histogram.
#' @param density Boolean dictating whether to slow density plot overlaid on histogram.
#' @param saveplot Boolean dictating whether to automatically save plot.
#' @param location Output location (file path) of saved plot (if any).
#' @return Plot of GP values.
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
      ggplot2::ggplot(data = gpvalues, ggplot2::aes(x = gp, y = ggplot2::after_stat(density))) + p +
        ggplot2::geom_density(alpha = 0, fill = "#000000", linewidth = 1) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = "GP Value", y = ylab) +
        ggplot2::ggtitle(titl) +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                       axis.title.x = ggplot2::element_text(size = 24), axis.title.y = ggplot2::element_text(size = 24),
                       axis.text.x = ggplot2::element_text(size = 20), axis.text.y = ggplot2::element_text(size = 20),
                       plot.title = ggplot2::element_text(size = 24, hjust = 0.5))
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

#' Ripley's H Function.
#'
#' Calculates the Ripley's H function.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y.
#' @param window Bounds of ROI to use, should be written as a vector: (xmin, xmax, ymin, ymax).
#' @param createplot Boolean dictating whether to create a plot of the Ripley's H function.
#' @param saveplot Boolean dictating whether to automatically save plot.
#' @param location Output location (file path) of saved plot (if any).
#' @return A matrix with two columns: the radii and the corresponding Ripley's H values.
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
            ggplot2::geom_line(ggplot2::aes(y = H), linewidth = 1) +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
            ggplot2::labs(x = "r", y = "H(r)") +
            ggplot2::ggtitle(" ") +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 24, hjust = 0.5),
                           panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                           axis.line = ggplot2::element_line(colour = "black"), axis.title.x = ggplot2::element_text(size = 32),
                           axis.text.x = ggplot2::element_text(size = 28), axis.text.y = ggplot2::element_text(size = 28),
                           axis.title.y = ggplot2::element_text(size = 32)))
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
#' Estimate the average cluster radius of a point cloud using the Ripley's H function.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y.
#' @return An estimate of the cluster radius.
#' @export
informradius <- function(points){
  #Get H values.
  Hr <- ripleyh(points, window = NULL, createplot = FALSE)
  #Find and return argmax.
  return(Hr[which.max(Hr[,2]),1])
}

#' Calculate intersections function used in Inform Persistence
#'
#' Calculates the point of intersection between two Gaussians.
#'
#' @param m1 Mean (peak) of first Gaussian.
#' @param m2 Mean (peak) of second Gaussian.
#' @param s1 Standard deviation of first Gaussian.
#' @param s2 Standard deviation of second Gaussian.
#' @param c1 Coefficient of first Gaussian.
#' @param c2 Coefficient of second Gaussian.
#' @return x value of point of intersection.
#' @export
calcinter <- function(m1, m2, s1, s2, c1, c2){
  a <- 1/(2 * s2 ** 2) - 1/(2 * s1 ** 2)
  b <- m1/(s1 ** 2) - m2/(s2 ** 2)
  c <- log(c1/c2) + log(s2/s1) + (m2 ** 2)/(2 * s2 ** 2) - (m1 ** 2)/(2 * s1 ** 2) 
  return(pracma::roots(c(a, b, c)))
}

#' Inform Deviance.
#'
#' Estimates the deviance parameter for JOSEPH clustering by fitting Gaussians to the distribution of marked values.
#'
#' @param values Vector of marked values from a marked point pattern (e.g. GP values).
#' @param k Number of Gaussians to fit. Should be equal to number of suspected categories.
#' @param createplot Boolean dictating whether to create a plot of the Gaussian fit.
#' @param saveplot Boolean dictating whether to automatically save plot.
#' @param location Output location (file path) of saved plot (if any).
#' @param outputall Boolean dictating whether to output all results of Gaussian fitting. Mostly used for development.
#' @param plotresolution Resolution of x axis used in plot (number of lines).
#' @return Output
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
          ggplot2::geom_line(ggplot2::aes(color = var), linewidth = 2) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::labs(x = "GP Value", y = "Density") +
          ggplot2::ggtitle("Estimated Distributions of GP Values") +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 24, hjust = 0.5),
                         panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(colour = "black"), axis.title.x = ggplot2::element_text(size = 24),
                         axis.text.x = ggplot2::element_text(size = 20), axis.text.y = ggplot2::element_text(size = 20),
                         axis.title.y = ggplot2::element_text(size = 24), legend.position = "none")
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

#Justification of Separation by Employed Persistent Homology

#' Justification of Separation by Employed Persistent Homology (JOSEPH)
#'
#' A topological data analysis tool for clustering marked point pattern data.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y, and the final column corresponding to mark.
#' @param deviance Deviance parameter - acceptable deviation from mean mark value of cluster.
#' @param radius Cluster radius.
#' @return The input marked point pattern with an additional column representing cluster allocation.
#' @export
joseph <- function(points, deviance = NULL, radius = NULL){
  #Check if points has too many columns.
  points <- points[,c(1,2,ncol(points))]
  #If undefined, calculate radius and deviance.
  if(is.null(radius)){
    radius <- informradius(points)
  }
  if(is.null(deviance)){
    deviance <- invisible(informdeviance(points[,3], createplot = FALSE))
    if(is.null(deviance)){
      #Standard deviation is probably 0.
      deviance <- min(abs(mean(points[,3]) - unique(points[,3])))
    }
  }
  
  #Calculate distance matrix.
  distance_matrix <- Rfast::Dist(points[,1:2])
  
  #Define number of points.
  number_of_points <- nrow(points)
  
  #Create repeating matrix of labels.
  label_matrix <- replicate(number_of_points, points[,3])
  #Exclude labels outside of radius for each point.
  distance_radius <- distance_matrix <= radius
  label_matrix[distance_matrix * distance_radius == 0] <- NA
  
  #Determine difference between density (mean value of all points in radius r) and value for each point.
  differences <- abs(points[,3] - colMeans(label_matrix, na.rm = TRUE))
  
  #Organise points in ascending order of difference.
  coords <- 1:number_of_points
  new_points <- as.data.frame(cbind(coords, differences))
  new_points <- new_points[order(new_points$differences),]
  
  #Extract the ordered differences. These will act as the coordinates.
  coordinates <- new_points[,1]
  
  #Extract filtrations.
  new_points <- cbind(new_points, 1:number_of_points)
  new_points <- new_points[order(new_points$coords),]
  filtrations <- new_points[,3]
  
  #Perform persistent homology.
  r <- rep(0, number_of_points)
  e <- rep(0, number_of_points)
  U <- list()
  
  #Iterate over all points starting from the lowest difference. Perform persistent homology to generate the partition.
  set_number <- 0
  for(i in coordinates){
    #Extract filtrations of neighbours of i.
    neighbour_filtrations <- filtrations[distance_radius[i,]]
    
    #Extract neighbours with lower filtration values than i.
    N <- coordinates[neighbour_filtrations[neighbour_filtrations < filtrations[i]]]
    
    #Check length of N.
    if(length(N) == 0){
      #If N is empty, then i is a peak. Create a new subset in U containing just i. Update r and e.
      r[i] <- i
      set_number <- set_number + 1
      e[i] <- set_number
      U[[set_number]] <- c(i)
    }else{
      #Organise the neighbours in N in descending order of the density of their roots.
      R <- r[N]
      neighbours <- as.data.frame(cbind(R, N))
      neighbours <- neighbours[order(neighbours$R),]
      
      #If N is non-empty, there must be at least one neighbour of higher density. Cycle through all neighbours.
      neighbour_differences <- abs(points[i,3] - points[neighbours[,1], 3])
      potentials <- neighbours[(1:nrow(neighbours))[neighbour_differences <= deviance], 1]
      if(length(potentials) == 0){
        #Node does not connect to any neighbours. Create new set for node.
        r[i] <- i
        set_number <- set_number + 1
        e[i] <- set_number
        U[[set_number]] <- c(i)
      } else {
        #Refine potentials.
        potentials <- unique(potentials)
        j <- potentials[1]
        U[[e[j]]] <- append(U[[e[j]]], i)
        r[i] <- j
        e[i] <- e[j]
        if(length(potentials) > 1){
          #Merge the sets of the potentials.
          for(k in potentials[-1]){
            e_k <- U[[e[k]]]
            U[[e[j]]] <- append(U[[e[j]]], U[[e[k]]])
            for(l in e_k){
              e[l] <- e[j]
              r[l] <- r[j]
            }
          }
        }
      }
    }
  }
  
  #Determine cluster indices.
  found <- 0
  #Filter out clusters with too few points.
  cluster_index <- rep(0, number_of_points)
  for(i in unique(r)){
    #Initialise new cluster.
    found <- found + 1
    cluster_index[U[[e[i]]]] <- found
  }
  #Return new point cloud with cluster index.
  return(cbind(points, cluster_index))
}

#' Filter JOSEPH.
#'
#' Filters JOSEPH cluster results based on user requirements.
#'
#' @param josephresults Output from JOSEPH clustering.
#' @param mode A string representing filter mode. Choose from "dispersion" (filter clusters with above-average dispersion), "minpoints" (removes clusters with point count lower than minpoints), "maxpoints" (removes clusters with point count higher than maxpoints) and "range" (removes clusters with mean marked value outside the range).
#' @param minpoints Minimum number of points.
#' @param maxpoints Maximum number of points.
#' @param range Vector representing lower and upper bounds of average mark to include.
#' @param radius Radius of clusters used to check dispersion.
#' @return Filtered JOSEPH results. The input marked point pattern with updated cluster allocation.
#' @export
filterjoseph <- function(josephresults, mode = "dispersion", minpoints = NULL, maxpoints = NULL, range = NULL, radius = NULL){
  #Determine number of clusters.
  #Check each possible mode.
  if("minpoints" %in% mode){
    #Check minpoints given.
    if(is.null(minpoints)){
      stop("No minimum given. Add argument minpoints = n, where n is the desired minimum number of points per cluster.")
    } else {
      #Get cluster numbers.
      cluster_numbers <- unique(josephresults[,4])
      cluster_numbers <- cluster_numbers[cluster_numbers > 0]
      #Iterate over all clusters.
      for(i in cluster_numbers){
        #Find points corresponding to that cluster.
        cluster_indicator <- josephresults[,4] == i
        if(sum(cluster_indicator) < minpoints){
          josephresults[cluster_indicator, 4] <- 0
        }
      }
    }
  }
  if("maxpoints" %in% mode){
    #Check minpoints given.
    if(is.null(maxpoints)){
      stop("No maximum given. Add argument maxpoints = n, where n is the desired maximum number of points per cluster.")
    } else {
      #Get cluster numbers.
      cluster_numbers <- unique(josephresults[,4])
      cluster_numbers <- cluster_numbers[cluster_numbers > 0]
      #Iterate over all clusters.
      for(i in cluster_numbers){
        #Find points corresponding to that cluster.
        cluster_indicator <- josephresults[,4] == i
        if(sum(cluster_indicator) > maxpoints){
          josephresults[cluster_indicator, 4] <- 0
        }
      }
    }
  }
  if("range" %in% mode){
    #Check range given.
    if(is.null(range)){
      stop("No range given. Add argument range = (m, n), where m and n are the lower and upper bounds of the mean labels per cluster, respectively.")
    } else {
      #Get cluster numbers.
      cluster_numbers <- unique(josephresults[,4])
      cluster_numbers <- cluster_numbers[cluster_numbers > 0]
      #Iterate over all clusters.
      for(i in cluster_numbers){
        #Find points corresponding to that cluster.
        cluster_indicator <- josephresults[,4] == i
        #Find mean value of that cluster.
        cluster_mean <- mean(josephresults[cluster_indicator,3])
        if(cluster_mean < range[1] | cluster_mean > range[2]){
          josephresults[cluster_indicator, 4] <- 0
        }
      }
    }
  }
  if("dispersion" %in% mode){
    #Normalise.
    original <- josephresults[,1:2]
    if(max(josephresults[,1]) - min(josephresults[,1]) > max(josephresults[,2]) - min(josephresults[,2])){
      josephresults[,1:2] <- (josephresults[,1:2] - min(josephresults[,1]))/(max(josephresults[,1]) - min(josephresults[,1]))
    } else {
      josephresults[,1:2] <- (josephresults[,1:2] - min(josephresults[,2]))/(max(josephresults[,2]) - min(josephresults[,2]))
    }
    #Calculate radius.
    if(is.null(radius)){
      #Calculate radius.
      radius <- informradius(josephresults)
    }
    #Calculate H(r) value.
    A <- (max(josephresults[,1]) - min(josephresults[,1])) * (max(josephresults[,2]) - min(josephresults[,2]))
    n <- nrow(josephresults)
    Hr <- sqrt(A) * sqrt((sum(Rfast::Dist(josephresults[,1:2]) <= radius) - n) / (2 * pi * n ** 2)) - radius
    #Get cluster numbers.
    cluster_numbers <- unique(josephresults[,4])
    cluster_numbers <- cluster_numbers[cluster_numbers > 0]
    #Iterate over all clusters.
    for(i in cluster_numbers){
      #Find points corresponding to that cluster.
      cluster_indicator <- josephresults[,4] == i
      #Extract points corresponding to that cluster.
      cluster_points <- josephresults[cluster_indicator, 1:2]
      #If only contains one point, delete as this will never show a large enough H value.
      if(nrow(cluster_points) <= 4){
        josephresults[cluster_indicator, 4] <- 0
      } else {
        #Calculate Ripley's H value at that radius for that cluster. If less than 0, remove cluster.
        n <- nrow(cluster_points)
        if(sqrt(A) * sqrt((sum(Rfast::Dist(cluster_points) <= radius) - n) / (2 * pi * n ** 2)) - radius < Hr){
          josephresults[cluster_indicator, 4] <- 0
        }
      }
    }
    #Return to original coordinates.
    josephresults[,1:2] <- original
  }
  #Return finished points.
  return(josephresults)
}

#' Split JOSEPH.
#'
#' Discretises JOSEPH results as in discretise(), with each cluster assigned a different category.
#'
#' @param josephresults Output from JOSEPH clustering.
#' @param splits Split values which define intervals to categorise points. If left as null, these intervals will be automatically estimated from Gaussian fitting.
#' @param k Number of Gaussians to fit in case splits is left null. Should be equal to number of suspected categories.
#' @param createplot Boolean dictating whether to create a plot of the new discretised marked point pattern.
#' @param location Output location (file path) of saved plot (if any).
#' @param saveplot Boolean dictating whether to automatically save plot.
#' @return A discretised marked point pattern.
#' @export
splitjoseph <- function(josephresults, splits = NULL, k = 2, createplot = FALSE, location = "", saveplot = FALSE){
  #Get original column names.
  originalnames <- colnames(josephresults)
  #Get splits if not specified.
  if(is.null(splits)){
    splits <- informdeviance(josephresults[,3], k = k, createplot = FALSE, outputall = TRUE)[[5]]
  }
  #Add extra column for categories.
  josephresults <- cbind(josephresults, rep(0, nrow(josephresults)))
  #Add infinities to split.
  splits <- c(-Inf, splits, Inf)
  #Determine number of clusters.
  cluster_numbers <- unique(josephresults[,4])
  #Iterate over each cluster.
  for(i in cluster_numbers[cluster_numbers != 0]){
    #Find points corresponding to that cluster.
    cluster_indicator <- josephresults[,4] == i
    #Extract cluster points.
    cluster_points <- josephresults[cluster_indicator,3]
    #Check number of points and calculate mean accordingly.
    if(length(cluster_points) == 1){
      cluster_mean <- cluster_points[1]
    } else {
      cluster_mean <- as.numeric(mean(cluster_points))
    }
    #Determine which category mean falls in.
    for(j in 1:(length(splits) - 1)){
      if(cluster_mean > splits[j] & cluster_mean <= splits[j + 1]){
        break
      }
    }
    #Update category column.
    josephresults[cluster_indicator,5] <- j
  }
  #Set background points to background.
  josephresults[josephresults[,4] == 0, 5] <- 0
  #Update column names.
  colnames(josephresults) <- c("x", "y", originalnames[3:4], "category")
  #Plot.
  if(createplot){
    print(
      p <- ggplot2::ggplot() +
        ggplot2::geom_point(data = josephresults, ggplot2::aes(x = x, y = y, colour = as.factor(category)), size = 2) +
        ggplot2::scale_color_brewer(palette = "Set2") +
        ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
        ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4))) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"),
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                       legend.position = "none") 
    )
    #Save plot.
    if(saveplot){
      ggplot2::ggsave(paste(location, "Split Results.pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
    }
  }
  #Return data.
  return(josephresults)
}

#' JOSEPH Cluster Summary function used in PLASMA.
#'
#' Outputs cluster statistics from JOSEPH results.
#'
#' @param josephresults Output from JOSEPH clustering.
#' @return A range of cluster statistics, including number of points, area of convex hull, density, average mark and mark standard deviation.
#' @export
josephclustersummary <- function(josephresults){
  #Check if string.
  if(is.character(josephresults)){
    return(rep(NULL, 7))
  } else {
    #Normalise.
    if(max(josephresults[,1]) - min(josephresults[,1]) > max(josephresults[,2]) - min(josephresults[,2])){
      scale <- max(josephresults[,1]) - min(josephresults[,1])
      josephresults[,1:2] <- (josephresults[,1:2] - min(josephresults[,1]))/(max(josephresults[,1]) - min(josephresults[,1]))
    } else {
      scale <- max(josephresults[,2]) - min(josephresults[,2])
      josephresults[,1:2] <- (josephresults[,1:2] - min(josephresults[,2]))/(max(josephresults[,2]) - min(josephresults[,2]))
    }
    #Get mean mark.
    meanmark <- mean(josephresults[,3])
    #Initialise empty data frame.
    clusterstatistics <- data.frame()
    #Calculate cluster numbers.
    cluster_numbers <- unique(josephresults[,4])
    if(length(cluster_numbers) == 1){
      print("JOSEPH: No clusters found.")
      return(rep(NULL, 7))
    }
    #Get number of clusters.
    cluster_numbers <- cluster_numbers[cluster_numbers != 0]
    #Iterate over each cluster number.
    for(i in cluster_numbers){
      #Extract points corresponding to that cluster.
      cluster_points <- josephresults[josephresults[,4] == i,]
      #Get area.
      pts <- cluster_points[chull(cluster_points[,1:2]),1:2]
      ar <- abs(pracma::polyarea(pts[,1], pts[,2])) * scale ** 2
      #Add statistics to data frame.
      clusterstatistics <- rbind(clusterstatistics, c(nrow(cluster_points), ar, nrow(cluster_points)/ar, mean(cluster_points[,3]), sd(cluster_points[,3]), meanmark, mean(cluster_points[,3]) - meanmark))
    }
    #Update column names.
    colnames(clusterstatistics) <- c("number of points", "area", "density", "average mark", "mark standard deviation", "global average", "difference")
    #Return cluster statistics.
    return(clusterstatistics)
  }
}

#' Plot Hulls.
#'
#' Plots convex hulls of clusters found by JOSEPH, coloured by mean GP value.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y, and the final column corresponding to mark.
#' @param josephresults Output from JOSEPH clustering.
#' @param mingp Minimum GP value for colour bar. Leave NULL for min in points.
#' @param maxgp Maximum GP value for colour bar. Leave NULL for max in points.
#' @param alpha Opacity alpha of hulls.
#' @param viridis Boolean dictating whether to use viridis colourscale - a cyan to magenta variant will be used if false.
#' @param saveplot Boolean dictating whether to automatically save plot.
#' @param location Output location (file path) of saved plot (if any).
#' @return Plot of hulls.
#' @export
plothulls <- function(points, josephresults, mingp = NULL, maxgp = NULL, alpha = 0.5, viridis = FALSE, saveplot = FALSE, location = ""){
  #Get min and max GP.
  if(is.null(mingp)){
    mingp <- min(points[,3])
  }
  if(is.null(maxgp)){
    maxgp <- max(points[,3])
  }
  means <- c()
  
  #Get min and max point coordinates
  minx <- min(points[,1])
  maxx <- max(points[,1])
  miny <- min(points[,2])
  maxy <- max(points[,2])
  
  #Create augmented data frame.
  data <- josephresults
  data <- data[data[,4] > 0,]
  
  #Initialise vertices data frame
  josephvers <- data.frame()
  count <- 0
  if(length(unique(data[,4])) > 0){
    for(i in unique(data[,4])){
      if(i > 0){
        count <- count + 1
        joe <- data[data[,4] == i,1:2]
        josephvers <- rbind(josephvers, cbind(joe[chull(joe),], count))
        means <- append(means, mean(data[data[,4] == i, 3]))
      }
    }
    
    #Get colours based on mean GPs.
    scaledmeans <- floor((means - mingp) / (maxgp - mingp) * 1000)
    if(viridis){
      vircols <- viridis::viridis(1000)
      cols <- vircols[scaledmeans]
    } else {
      colramp <- colorRampPalette(c("cyan", "magenta"))
      cols <- colramp(1000)[scaledmeans]
    }
    
    
    #Create plot.
    colnames(josephvers) <- c("x", "y", "v")
    print(ggplot2::ggplot(josephvers, ggplot2::aes(x = x, y = y, fill = factor(v))) +
            ggplot2::geom_polygon(alpha = alpha) +
            ggplot2::scale_x_continuous(limits = c(minx, maxx), expand = c(0.01, 0.01)) +
            ggplot2::scale_y_continuous(limits = c(miny, maxy), expand = c(0.01, 0.01)) +
            ggplot2::scale_fill_manual(values = cols) +
            ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                           axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                           legend.position = "none"))
  }
  #Save plot.
  if(saveplot){
    ggplot2::ggsave(paste(location, "Hulls.svg", sep = ""), device = "svg", width = 14.342, height = 14.342, units = "cm")
  }
}

#' Compute Intersection Over Union.
#'
#' Takes the labelled vertices of two sets of convex hulls. Calculates the area of intersection over the area of the union of the hulls.
#'
#' @param vertices_1 Vertices of first set of hulls. Represented as a matrix with three columns: x, y and vertex label.
#' @param vertices_2 Vertices of second set of hulls. Represented as a matrix with three columns: x, y and vertex label.
#' @return IOU score.
#' @export
computeIOU <- function(vertices_1, vertices_2){
  #Convert to data frames.
  vertices_1 <- as.data.frame(vertices_1)
  vertices_2 <- as.data.frame(vertices_2)
  
  #Fix column names.
  colnames(vertices_1) <- c("x", "y", "hull")
  colnames(vertices_2) <- c("x", "y", "hull")
  #Calculate convex hulls and convert to polygons.
  hulls_1 <- vertices_1 %>% 
    dplyr::group_by(hull) %>% 
    dplyr::summarise(geometry = list(sf::st_polygon(list(cbind(x[c(chull(x, y), chull(x, y)[1])], y[c(chull(x, y), chull(x, y)[1])]))))) %>%
    dplyr::summarise(geometry = list(do.call(sf::st_union, geometry)))
  
  hulls_2 <- vertices_2 %>% 
    dplyr::group_by(hull) %>% 
    dplyr::summarise(geometry = list(sf::st_polygon(list(cbind(x[c(chull(x, y), chull(x, y)[1])], y[c(chull(x, y), chull(x, y)[1])]))))) %>%
    dplyr::summarise(geometry = list(do.call(sf::st_union, geometry)))
  
  #Convert to sf objects.
  sf_hulls_1 <- sf::st_as_sf(hulls_1)
  sf_hulls_2 <- sf::st_as_sf(hulls_2)
  
  #Calculate intersection and union of all hulls.
  intersection <- sf::st_intersection(sf_hulls_1$geometry[[1]], sf_hulls_2$geometry[[1]])
  union <- sf::st_union(sf_hulls_1$geometry[[1]], sf_hulls_2$geometry[[1]])
  
  #Calculate IoU
  iou <- ifelse(sf::st_is_empty(intersection), 0, sf::st_area(intersection) / sf::st_area(union))
  return(iou)
}

#' Compute IOU from JOSEPH.
#'
#' Computes the IOU between a defined set of hulls and the hulls identified from JOSEPH clustering.
#'
#' @param vertices_1 Vertices of first set of hulls. Represented as a matrix with three columns: x, y and vertex label.
#' @param josephresults Output from JOSEPH clustering.
#' @return IOU score.
#' @export
computeIOUfromjoseph <- function(vertices_1, josephresults){
  #Convert to data frames.
  vertices_1 <- as.data.frame(vertices_1)
  josephresults <- as.data.frame(josephresults)
  
  #Fix column names.
  colnames(vertices_1) <- c("x", "y", "hull")
  
  #Calculate cluster hull vertices.
  res <- josephresults[josephresults[,4] > 0, c(1,2,4)]
  vertices_2 <- do.call(rbind, lapply(unique(res[,3]), function(i){
    joe <- res[res[,3] == i, 1:2]
    return(cbind(joe[chull(joe), 1:2], i))
  }))
  colnames(vertices_2) <- c("x", "y", "hull")
  
  #Calculate convex hulls and convert to polygons.
  hulls_1 <- vertices_1 %>% 
    dplyr::group_by(hull) %>% 
    dplyr::summarise(geometry = list(sf::st_polygon(list(cbind(x[c(chull(x, y), chull(x, y)[1])], y[c(chull(x, y), chull(x, y)[1])]))))) %>%
    dplyr::summarise(geometry = list(do.call(sf::st_union, geometry)))
  
  hulls_2 <- vertices_2 %>% 
    dplyr::group_by(hull) %>% 
    dplyr::summarise(geometry = list(sf::st_polygon(list(cbind(x[c(chull(x, y), chull(x, y)[1])], y[c(chull(x, y), chull(x, y)[1])]))))) %>%
    dplyr::summarise(geometry = list(do.call(sf::st_union, geometry)))
  
  #Convert to sf objects.
  sf_hulls_1 <- sf::st_as_sf(hulls_1)
  sf_hulls_2 <- sf::st_as_sf(hulls_2)
  
  #Calculate intersection and union of all hulls.
  intersection <- sf::st_intersection(sf_hulls_1$geometry[[1]], sf_hulls_2$geometry[[1]])
  union <- sf::st_union(sf_hulls_1$geometry[[1]], sf_hulls_2$geometry[[1]])
  
  #Calculate IoU
  iou <- ifelse(sf::st_is_empty(intersection), 0, sf::st_area(intersection) / sf::st_area(union))
  return(iou)
}

#' Compute IOU from JOSEPHs.
#'
#' Computes the IOU between two sets of hulls identified from JOSEPH clustering.
#'
#' @param josephresults_1 Output from first JOSEPH clustering.
#' @param josephresults_2 Output from second JOSEPH clustering.
#' @return IOU score.
#' @export
computeIOUfromjosephs <- function(josephresults_1, josephresults_2){
  #Convert to data frames.
  josephresults_1 <- as.data.frame(josephresults_1)
  josephresults_2 <- as.data.frame(josephresults_2)
  
  #Check for clusters.
  if(max(josephresults_1[,4]) == 0 || max(josephresults_2[,4]) == 0){
    return(0)
  }
  
  #Calculate cluster hull vertices.
  res <- josephresults_1[josephresults_1[,4] > 0, c(1,2,4)]
  vertices_1 <- do.call(rbind, lapply(unique(res[,3]), function(i){
    joe <- res[res[,3] == i, 1:2]
    return(cbind(joe[chull(joe), 1:2], i))
  }))
  colnames(vertices_1) <- c("x", "y", "hull")
  
  #Calculate cluster hull vertices.
  res <- josephresults_2[josephresults_2[,4] > 0, c(1,2,4)]
  vertices_2 <- do.call(rbind, lapply(unique(res[,3]), function(i){
    joe <- res[res[,3] == i, 1:2]
    return(cbind(joe[chull(joe), 1:2], i))
  }))
  colnames(vertices_2) <- c("x", "y", "hull")
  
  #Calculate convex hulls and convert to polygons.
  hulls_1 <- vertices_1 %>% 
    dplyr::group_by(hull) %>% 
    dplyr::summarise(geometry = list(sf::st_polygon(list(cbind(x[c(chull(x, y), chull(x, y)[1])], y[c(chull(x, y), chull(x, y)[1])]))))) %>%
    dplyr::summarise(geometry = list(do.call(sf::st_union, geometry)))
  
  hulls_2 <- vertices_2 %>% 
    dplyr::group_by(hull) %>% 
    dplyr::summarise(geometry = list(sf::st_polygon(list(cbind(x[c(chull(x, y), chull(x, y)[1])], y[c(chull(x, y), chull(x, y)[1])]))))) %>%
    dplyr::summarise(geometry = list(do.call(sf::st_union, geometry)))
  
  #Convert to sf objects.
  sf_hulls_1 <- sf::st_as_sf(hulls_1)
  sf_hulls_2 <- sf::st_as_sf(hulls_2)
  
  #Calculate intersection and union of all hulls.
  intersection <- sf::st_intersection(sf_hulls_1$geometry[[1]], sf_hulls_2$geometry[[1]])
  union <- sf::st_union(sf_hulls_1$geometry[[1]], sf_hulls_2$geometry[[1]])
  
  #Calculate IoU
  iou <- ifelse(sf::st_is_empty(intersection), 0, sf::st_area(intersection) / sf::st_area(union))
  return(iou)
}

#' Simple Hull Plot.
#'
#' Plots the hulls from two sets of hulls.
#'
#' @param vertices_1 Vertices of first set of hulls. Represented as a matrix with three columns: x, y and vertex label.
#' @param vertices_2 Vertices of second set of hulls. Represented as a matrix with three columns: x, y and vertex label.
#' @param legend Boolean dictating whether to show legend in plot.
#' @param legend_title String which represents the title of the legend (if shown).
#' @param legend_labels Vector of strings representing labels of hull types (leave NULL for numeric default).
#' @param limits_x Vector representing the minimum and maximum x value of the ROI.
#' @param limits_y Vector representing the minimum and maximum y value of the ROI.
#' @return Output
#' @export
simplehullplot <- function(vertices_1, vertices_2, legend = FALSE, legend_title = " ", legend_labels = NULL, limits_x = NULL, limits_y = NULL){
  #Add group column.
  if(nrow(vertices_1) > 0){
    vertices_1 <- cbind(vertices_1, 1)
    colnames(vertices_1) <- c("x", "y", "hull", "group")
    vertmax <- max(vertices_1$hull)
  } else {
    vertmax <- 0
    no_1 <- TRUE
  }
  if(nrow(vertices_2) > 0){
    vertices_2 <- cbind(vertices_2, 2)
    colnames(vertices_2) <- c("x", "y", "hull", "group")
  } else {
    no_2 <- TRUE
  }
  
  #Check if data is available.
  if(no_1 && no_2){
    return(FALSE)
  }
  
  #Combine the hulls into one list
  vertices_2$hull <- vertices_2$hull + vertmax
  
  #Calculate the convex hulls
  hulls_1 <- lapply(split(vertices_1, vertices_1$hull), function(df) df[chull(df$x, df$y), ])
  hulls_2 <- lapply(split(vertices_2, vertices_2$hull), function(df) df[chull(df$x, df$y), ])
  hulls_df <- rbind(do.call(rbind, hulls_1), do.call(rbind, hulls_2))
  
  #Create plot.
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = hulls_df, ggplot2::aes(x = x, y = y, fill = as.factor(group), group = hull), alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "X", y = "Y", fill = "Group") +
    ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                   axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())
  
  #Set boundaries.
  if(is.null(limits_x)){
    p <- p + ggplot2::scale_x_continuous(expand = c(0.01, 0.01))
  } else {
    p <- p + ggplot2::scale_x_continuous(limits = limits_x, expand = c(0.01, 0.01))
  }
  if(is.null(limits_y)){
    p <- p + ggplot2::scale_y_continuous(expand = c(0.01, 0.01))
  } else {
    p <- p + ggplot2::scale_y_continuous(limits = limits_y, expand = c(0.01, 0.01))
  }
  
  #Legend.
  if(legend && !is.null(legend_labels)){
    p <- p + ggplot2::scale_fill_discrete(name = legend_title, labels = legend_labels)
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  #Print plot.
  print(p)
  return(TRUE)
}

#Marked Point Pattern Histogram.

#' Bin Datum function used in MPP Histogram
#'
#' Bins data points together. Dev use only.
#'
#' @param i Bin number.
#' @param y_wid Width of bin in y direction.
#' @param GP_x GP in x.
#' @param GP_y GP in y.
#' @param pts Points.
#' @return Binned data.
#' @export
bindatum <- function(i, y_wid, GP_x, GP_y, pts){
  return(c(1 + (i - 1) %/% y_wid, 1 + (i - 1) %% y_wid,
           sum(pts[GP_x == 1 + (i - 1) %/% y_wid & GP_y == 1 + (i - 1) %% y_wid, 3]),
           sum(pts[GP_x == 1 + (i - 1) %/% y_wid & GP_y == 1 + (i - 1) %% y_wid, 4])))
}

#' Compile Matrix function used in MPP Histogram.
#'
#' Converts bindatum output into matrix form. Dev use only. 
#'
#' @param i Bin number.
#' @param rep Representative.
#' @param col Column number.
#' @return Matrix form of bin data.
#' @export
compilematrix <- function(i, rep, col){
  return(rep[rep[,1] == i, col])
}

#' MPP Histogram.
#'
#' Creates Membrane Order Map (MOM).
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y, and the final column corresponding to mark.
#' @param Width of bins to use in histogram.
#' @param location Save location (file path) of plots.
#' @param plotname Name of plot for file name.
#' @param useggplot Boolean dictating whether to use ggplot to construct histogram.
#' @param outputdata Boolean dictating whether to output all results (dev use only).
#' @return A membrane order map.
#' @export
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

#Marked Point Pattern Simulator.

#' Toroidal Distance Matrix function used in MPP Simulate.
#'
#' Calculates toroidal distance matrix of a point cloud.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y.
#' @param domain_max_x Size of ROI in x direction.
#' @param domain_max_y Size of ROI in y direction.
#' @return Toroidal distance matrix.
#' @export
toroidaldist <- function(points, domain_max_x, domain_max_y){
  #Compute x and y differences.
  x_d <- do.call(cbind, lapply(1:nrow(points), function(i){
    x_d <- abs(points[i,1] - points[,1])
    x_d[x_d > domain_max_x / 2] <-  domain_max_x - x_d[x_d > domain_max_x / 2]
    return(x_d)
  }))
  y_d <- do.call(cbind, lapply(1:nrow(points), function(i){
    y_d <- abs(points[i,2] - points[,2])
    y_d[y_d > domain_max_y / 2] <-  domain_max_y - y_d[y_d > domain_max_y / 2]
    return(y_d)
  }))
  #Create toroidal coordinates matrix.
  return(sqrt(x_d ** 2 + y_d ** 2))
}

#' Calculate Getis L Values function used in MPP Simulate.
#'
#' Calculates Getis L values (dev use only).
#'
#' @param distance_matrix Distance matrix of point cloud.
#' @param number_of_points Number of points in point cloud.
#' @param domain_size_x Size of ROI in x direction.
#' @param domain_size_y Size of ROI in y direction.
#' @param search_radius Radius to calculate L over.
#' @return L value for each point.
#' @export
calculate_L_values <- function(distance_matrix, number_of_points, domain_size_x, domain_size_y, search_radius){
  return(do.call(c, lapply(1:number_of_points, function(i){
    return(sqrt(((domain_size_x * domain_size_y) / (pi * number_of_points)) * (sum(rep(1,number_of_points)[distance_matrix[i,] <= search_radius])-1)))
  })))
}

#' Velocity function used in MPP Simulate.
#'
#' Determines per-frame velocity of each point (dev use only).
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y.
#' @param number_of_points Number of points in point cloud.
#' @param domain_size_x Size of ROI in x direction.
#' @param domain_size_y Size of ROI in y direction.
#' @param search_radius Radius to calculate L over.
#' @param D_max Max step size.
#' @param L_target Target L value.
#' @return Point offsets.
#' @export
velocity <- function(points, number_of_points, domain_size_x, domain_size_y, search_radius, D_max, L_target){
  #Calculate L values.
  L_values <- calculate_L_values(distance_matrix = toroidaldist(points, domain_max_x = domain_size_x/2, domain_max_y = domain_size_y/2), number_of_points, domain_size_x, domain_size_y, search_radius)
  #Apply quadratically decreasing step size dependent on Ripley's value.
  return(do.call(rbind, lapply(1:number_of_points, function(i){
    angle <- runif(n = 1, min = 0, max = 2 * pi)
    step <- D_max * (1 - (L_values[i]/L_target) ** 2)
    if(step < 0){
      v_x <- 0
      v_y <- 0
    } else{
      v_x <- step * cos(angle)
      v_y <- step * sin(angle)
    }
    return(c(v_x, v_y))
  })))
}

#' MPP Simulate.
#'
#' Versatile simulator for generating marked point patterns with pre-defined properties. Produces spatial clusters by trying to fit to a given Getis L value over a given radius. These are then overlaid with marked point clusters. Each marked point is randomly assigned a category, with corresponding Gaussian distribution for generating marked values.
#'
#' @param number_of_frames Number of frames to simulate over.
#' @param number_of_points Number of points to use in simulation.
#' @param domain_size_x Size of ROI in x direction.
#' @param domain_size_y Size of ROI in y direction. Will be same as x if left NULL.
#' @param search_radius Radius to calculate Getis L value.
#' @param L_target Target L value to fit to.
#' @param D_max Max step size between frames.
#' @param tomato_threshold Persistence threshold for ToMATo clustering when producing spatial clusters. Leave NULL for automatic estimation.
#' @param grid_size Granularity of grid used to place marked points overlaid on spatial clusters.
#' @param number_of_categories Number of marked point types, each with own Normal distribution of marked values.
#' @param number_of_resimulated_points Number of points to use in marked point pattern, leave as NULL for same as number_of_points.
#' @param use_csr Boolean dictating whether to use a completely spatially random distribution for the marked point pattern.
#' @param probability_of_picking_categories Probability of assigning each spatial cluster to a given category.
#' @param precision_alpha Alpha parameter of the Gamma distribution used to simulate localisation error.
#' @param precision_beta Beta parameter of the Gamma distribution used to simulate localisation error.
#' @param category_mean Vector of mean values of the Normal distributions used to generate marks for each category.
#' @param category_sd Vector of standard deviations of the Normal distributions used to generate marks for each category.
#' @param background_mean Mean value of the Normal distribution used to generate marks for the background points.
#' @param background_sd Standard deviation of the Normal distribution used to generate marks for the background points.
#' @param createplots Boolean dictating whether to create plots of the marked point pattern and hull heatmap.
#' @param saveplots Boolean dictating whether to automatically save plots.
#' @param location Output location (file path) of saved plots (if any).
#' @param plotname Name of plot for file name.
#' @param heatmapname Name of heatmap for file name.
#' @param outputall Boolean dictating whether all results should be output (dev use only).
#' @return A simulated marked point pattern. A matrix with 4 columns: x, y, assigned category and marked value.
#' @export
mppsimulate <- function(number_of_frames = 100, number_of_points = 500, domain_size_x = 1000, domain_size_y = NULL,
                        search_radius = 50, L_target = 100, D_max = 50, tomato_threshold = NULL, grid_size = 100,
                        number_of_categories = 1, number_of_resimulated_points = NULL, use_csr = FALSE,
                        probability_of_picking_categories = c(0.7), precision_alpha = 0, precision_beta = 0,
                        category_mean = c(1), category_sd = c(0), background_mean = 0, background_sd = 0,
                        createplots = FALSE, saveplots = FALSE, location = "", plotname = "Plot",
                        heatmapname = "Heatmap", outputall = FALSE){
  #Get domain size.
  if(is.null(domain_size_y)){
    domain_size_y <- domain_size_x
  }
  
  #Define number of points to simulate.
  if(is.null(number_of_resimulated_points)){
    number_of_resimulated_points <- number_of_points
  }
  
  #Check radius.
  if(search_radius < 30){
    print("MPP Simulate: Use a larger search radius to produce noticeable clusters. Setting radius to 30.")
    search_radius <- 30
  }
  
  #Simulate uniformly distributed initial points
  points <- cbind(points_x = runif(n = number_of_points, min = 0, max = domain_size_x),
                  points_y = runif(n = number_of_points, min = 0, max = domain_size_y))
  
  #Repeat until all frames have been produced.
  for(i in 1:(number_of_frames - 1)){
    #Add velocity.
    points <- points + velocity(points, number_of_points, domain_size_x, domain_size_y, search_radius, D_max, L_target)
    #Return points to domain.
    points[points[,1] < 0, 1] <- points[points[,1] < 0, 1] + domain_size_x
    points[points[,1] > domain_size_x, 1] <- points[points[,1] > domain_size_x, 1] - domain_size_x
    points[points[,2] < 0, 2] <- points[points[,2] < 0, 2] + domain_size_y
    points[points[,2] > domain_size_y, 2] <- points[points[,2] > domain_size_y, 2] - domain_size_y
  }
  
  #Create domain field.
  #Domain field parameters.
  grid_length_x <- domain_size_x / grid_size
  grid_length_y <- domain_size_y / grid_size
  
  #Split field into grid
  field <- matrix(0L, nrow = grid_size, ncol = grid_size)
  field_data <- data.frame()
  
  #Run ToMATo on the simulated data to determine approximate coordinates
  #Set ToMATo parameters
  tomato_radius <- search_radius
  
  #If threshold not given, calculate.
  if(is.null(tomato_threshold)){
    tomato_threshold <- informclusterpersistence(points, radius = tomato_radius)
  }
  #If still NULL, input data showed no clustering.
  if(is.null(tomato_threshold)){
    print("MPP Simulate: The input data shows no sign of clustering, so we cannot simulate domains. Try adjusting the search radius, L target value or number of points.")
  }
  
  #Cluster data using ToMATo and determine number of clusters
  cluster_index <- getari(points, tomato_radius, tomato_threshold)[,3]
  number_of_clusters <- sum(unique(cluster_index) > 0)
  
  #For each cluster (if there are any), determine the convex hull and set the corresponding squares in the field
  #matrix to be domains.
  domain_vertices <- list()
  if(number_of_clusters > 0){
    #Create centre points data frame.
    centre_points <- as.data.frame(matrix(0, nrow = grid_size ** 2, ncol = 2))
    for(j in 1:grid_size){
      centre_points[((j - 1) * grid_size + 1):(j * grid_size), 1:2] <- cbind(rep((j - 0.5) * grid_length_x, grid_size), rep(j, grid_size))
    }
    centre_points <- cbind(centre_points, (1:grid_size - 0.5) * grid_length_y, 1:grid_size, 0)
    #Check if each centre point is inside each cluster.
    for(i in 1:number_of_clusters){
      #Extract cluster points.
      cluster_coordinates <- points[cluster_index == i, 1:2]
      #Filter out clusters with less than three points.
      if(length(cluster_coordinates) > 3){
        #Determine convex hull.
        hull_points <- chull(cluster_coordinates)
        vertices <- cbind(cluster_coordinates[hull_points, 1], cluster_coordinates[hull_points, 2])
        domain_vertices[[length(domain_vertices) + 1]] <- vertices
        #Check if the centre point of each square is within the convex hull.
        centre_points[pracma::inpolygon(centre_points[,1], centre_points[,3], vertices[,1], vertices[,2]),5] <- i
      }
    }
    #Create field data.
    field_data <- centre_points[,c(2,4,5)]
    #Filter out background points.
    centre_points <- centre_points[centre_points[,5] != 0,]
    #Update field.
    for(i in 1:nrow(centre_points)){
      field[centre_points[i,2], centre_points[i,4]] <- centre_points[i,5]
    }
    
    #Assign each square to a category, amend this to field data frame.
    categories <- 1:number_of_categories
    dif <- number_of_clusters - number_of_categories
    if(dif < 0){
      #If there are more categories than domains, assign as many domains to unique categories as possible.
      domain_categories <- append(0, sample(categories, number_of_clusters))
    } else{
      #Randomly label each domain with a certain category.
      domain_categories <- append(0, sample(c(categories, sample(categories, size = dif, replace = TRUE)), size = number_of_clusters))
    }
    
    #Record the category allocation of each square.
    field_data[,3] <- field_data[,3] + 1
    field_data <- cbind(field_data, 0, grid_size - field_data[,1] + 1)
    field_data[,4] <- domain_categories[field_data[,3]]
    
    #If using CSR, create probability vector
    if(use_csr){
      probability_of_picking_categories <- c()
      number_of_grids <- grid_size * grid_size
      for(i in categories){
        probability_of_picking_categories <- append(probability_of_picking_categories, length(rep(0, number_of_grids)[field_data[,4] == i])/number_of_grids)
      }
    }
    
    #Change column names in field data frame
    colnames(field_data) <- c("x","y","label", "cat", "invx")
    
    #Generate points sequentially.
    points_x <- rep(0, number_of_resimulated_points)
    points_y <- rep(0, number_of_resimulated_points)
    point_labels <- rep(0, number_of_resimulated_points)
    original_category <- rep(0, number_of_resimulated_points)
    all_categories <- append(0, categories)
    all_probabilities <- append(1 - sum(probability_of_picking_categories), probability_of_picking_categories)
    number_of_category_squares <- nrow(field_data)
    for(i in 1:number_of_resimulated_points){
      #Assign the point a category and record that category
      point_category <- sample(all_categories, size = 1, prob = all_probabilities)
      original_category[i] <- point_category
      #Pick a square in that category
      square_number <- sample((1:number_of_category_squares)[field_data[,4] == point_category], size = 1)
      #Generate the point at a random location in that square
      square_i <- field_data[square_number, 5]
      square_j <- field_data[square_number, 2]
      points_x[i] <- runif(n = 1, min = grid_length_x * (square_j - 1), max = grid_length_x * square_j)
      points_y[i] <- domain_size_y - runif(n = 1, min = grid_length_y * (square_i - 1), max = grid_length_y * square_i)
      #Assign the point a value
      if(point_category == 0){
        point_labels[i] <- rnorm(1, mean = background_mean, sd = background_sd)
      } else{
        point_labels[i] <- rnorm(1, mean = category_mean[point_category], sd = category_sd[point_category])
      }
    }
  } else {
    points_x <- runif(number_of_resimulated_points, min = 0, max = domain_size_x)
    points_y <- runif(number_of_resimulated_points, min = 0, max = domain_size_y)
    point_labels <- rnorm(number_of_resimulated_points, mean = background_mean, sd = background_sd)
    original_category <- rep(background_mean, number_of_resimulated_points)
    #Create empty field data frame.
    field_data <- c()
    for(i in 1:grid_size){
      field_data <- c(field_data, rep(i, grid_size))
    }
    field_data <- as.data.frame(cbind(field_data, grid_size:1, 1, 0))
    #Change column names in field data frame
    colnames(field_data) <- c("x","y","label","cat")
  }
  
  #Combine all points data.
  points <- as.data.frame(cbind(points_x, points_y, point_labels, original_category))
  
  #Apply localisation precision
  if(precision_alpha > 0 && precision_beta > 0){
    #Create and save plots original plots.
    if(createplots){
      #Check for discrete case and plot accordingly.
      if(length(category_sd) == length(category_sd[category_sd == 0])){
        #Plot labelled data.
        print(
          ggplot2::ggplot(data = points, ggplot2::aes(x = points_x, y = points_y)) +
            ggplot2::geom_point(ggplot2::aes(color = factor(point_labels))) +
            ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
            ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
            ggplot2::labs(color = "Value") +
            ggplot2::scale_colour_brewer(palette = "Set2") +
            ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                           axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                           legend.title = ggplot2::element_text(size = 20), legend.text = ggplot2::element_text(size = 16))
        )
        #Save plot.
        if(saveplots){
          ggplot2::ggsave(paste(location, plotname, " (Pre-localisation Error).pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
        }
      } else{
        #Plot labelled data in continuous case.
        print(
          ggplot2::ggplot(data = points, ggplot2::aes(x = points_x, y = points_y)) +
            ggplot2::geom_point(ggplot2::aes(color = point_labels)) +
            ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
            ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                           axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                           legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.75, "cm"),
                           legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16))
        )
        #Save plot.
        if(saveplots){
          ggplot2::ggsave(paste(location, plotname, " (Pre-localisation Error).pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
        }
      }
    }
    #Simulate sds from Gamma distribution.
    sds <- rgamma(n = number_of_resimulated_points, shape = precision_alpha, scale = precision_beta)
    #Simulate angles.
    angles <- runif(n = number_of_resimulated_points, min = 0, max = 2 * pi)
    #Create distances.
    distances <- rnorm(n = number_of_resimulated_points, mean = 0, sd = sds)
    #Update point positions.
    points[,1] <- points[,1] + distances * cos(angles)
    points[,2] <- points[,2] + distances * sin(angles)
    #Apply toroidal coordinate system.
    points[points[,1] < 0, 1] <- points[points[,1] < 0, 1] + domain_size_x
    points[points[,1] > domain_size_x, 1] <- points[points[,1] > domain_size_x, 1] - domain_size_x
    points[points[,2] < 0, 2] <- points[points[,2] < 0, 2] + domain_size_y
    points[points[,2] > domain_size_y, 2] <- points[points[,2] > domain_size_y, 2] - domain_size_y
  }
  
  #Create plots.
  if(createplots){
    #Check for discrete case and plot accordingly.
    if(length(category_sd) == length(category_sd[category_sd == 0])){
      #Plot labelled data.
      print(
        ggplot2::ggplot(data = points, ggplot2::aes(x = points_x, y = points_y)) +
          ggplot2::geom_point(ggplot2::aes(color = factor(point_labels))) +
          ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
          ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
          ggplot2::labs(color = "Value") +
          ggplot2::scale_colour_brewer(palette = "Set2") +
          ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                         axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                         legend.title = ggplot2::element_text(size = 20), legend.text = ggplot2::element_text(size = 16))
      )
      #Save plot.
      if(saveplots){
        ggplot2::ggsave(paste(location, plotname, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
      }
    } else{
      #Plot labelled data in continuous case.
      print(
        ggplot2::ggplot(data = points, ggplot2::aes(x = points_x, y = points_y)) +
          ggplot2::geom_point(ggplot2::aes(color = point_labels)) +
          ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
          ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                         axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                         legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.75, "cm"),
                         legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16))
      )
      #Save plot.
      if(saveplots){
        ggplot2::ggsave(paste(location, plotname, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
      }
    }
    
    #Create and save heatmap.
    print(
      ggplot2::ggplot(field_data) + ggplot2::aes(x = y, y = x, fill = as.factor(cat)) + ggplot2::geom_tile() +
        ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
        ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
        ggplot2::labs(fill = "Category") +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(), legend.title = ggplot2::element_text(size = 20),
                       legend.text = ggplot2::element_text(size = 16))
    )
    #Save heatmap.
    if(saveplots){
      ggplot2::ggsave(paste(location, heatmapname, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
    }
  }
  
  #Return points.
  if(outputall){
    if(length(category_mean) == 1){
      if(category_mean == background_mean){
        overlap <- Inf
      } else {
        overlap <- abs(sqrt(category_sd ** 2 + background_sd ** 2) / (category_mean - background_mean))
      }
      return(list(points = points, 
                  inputparameters = cbind(c("number_of_frames", "number_of_points", "number_of_resimulated_points", "domain_size_x", "domain_size_y", "search_radius", "L_target", "D_max", "tomato_threshold", "grid_size", "number_of_categories", "use_csr", "probability_of_picking_categories", "precision_alpha", "precision_beta", "category_mean", "category_sd", "background_mean", "background_sd"),
                                          c(number_of_frames, number_of_points, number_of_resimulated_points, domain_size_x, domain_size_y, search_radius, L_target, D_max, tomato_threshold, grid_size, number_of_categories, use_csr, probability_of_picking_categories, precision_alpha, precision_beta, category_mean, category_sd, background_mean, background_sd)),
                  summary = cbind(c("number_of_domains", "area_of_domains", "mean_localisation_error", "mark_overlap"),
                                  c(number_of_clusters, domain_size_x * domain_size_y * sum(field_data[,4] > 0) / nrow(field_data), precision_alpha * precision_beta, overlap)),
                  domainvertices = domain_vertices))
    } else {
      return(list(points = points, 
                  inputparameters = cbind(c("number_of_frames", "number_of_points", "number_of_resimulated_points", "domain_size_x", "domain_size_y", "search_radius", "L_target", "D_max", "tomato_threshold", "grid_size", "number_of_categories", "use_csr", "probability_of_picking_categories", "precision_alpha", "precision_beta", "category_mean", "category_sd", "background_mean", "background_sd"),
                                          c(number_of_frames, number_of_points, number_of_resimulated_points, domain_size_x, domain_size_y, search_radius, L_target, D_max, tomato_threshold, grid_size, number_of_categories, use_csr, paste(probability_of_picking_categories, collapse = " "), precision_alpha, precision_beta, paste(category_mean, collapse = " "), paste(category_sd, collapse = " "), background_mean, background_sd)),
                  summary = cbind(c("number_of_domains", "area_of_domains", "mean_localisation_error"),
                                  c(number_of_clusters, domain_size_x * domain_size_y * sum(field_data[,4] > 0) / nrow(field_data), precision_alpha * precision_beta)),
                  domainvertices = domain_vertices))
    }
  } else {
    return(points)
  }
}

#' Segment function used in MPP Segment.
#'
#' Segments an MPP (dev use only).
#'
#' @param x x coordinate.
#' @param interval Interval to segment.
#' @param coords Coordinates of points.
#' @return Segmented interval.
#' @export
seg <- function(x, interval, coords){
  #Take the given interval and extract only points in that interval.
  return(coords[coords[,1] >= interval[x,1] & coords[,1] < interval[x,2] &
                  coords[,2] >= interval[x,3] & coords[,2] < interval[x,4],])
}

#' Intervals function used in MPP Segment.
#'
#' Segments an MPP (dev use only).
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y, and the final column corresponding to mark.
#' @param m Number of x intervals.
#' @param n Number of y intervals. Will use m if left NULL.
#' @param xlims Vector of minimum and maximum x values to use.
#' @param ylims Vector of minimum and maximum y values to use.
#' @return Segmented intervals.
#' @export
intervals <- function(points, m, n = NULL, xlims = NULL, ylims = NULL){
  #Determine NULL parameters.
  if(is.null(n)){
    n <- m
  }
  #Impose limits.
  if(!is.null(xlims)){
    points <- points[points[,1] >= xlims[1] & points[,1] <= xlims[2],]
  } else {
    xlims <- c(min(points[,1]), max(points[,1]))
  }
  if(!is.null(ylims)){
    points <- points[points[,2] >= ylims[1] & points[,2] <= ylims[2],]
  } else {
    ylims <- c(min(points[,2]), max(points[,2]))
  }
  
  #Bring to origin.
  points[,1] <- points[,1] - xlims[1]
  points[,2] <- points[,2] - ylims[1]
  
  #Define segment size in x and y directions.
  xseg <- (xlims[2] - xlims[1])/m
  yseg <- (ylims[2] - ylims[1])/n
  
  #Create intervals.
  interv <- matrix(0L, nrow = m * n, ncol = 4)
  for(i in 1:m){
    interv[((i - 1) * n + 1):(i * n),] <- cbind(rep(i - 1, n) * xseg, rep(i, n) * xseg, 0:(n-1) * yseg, 1:n * yseg)
  }
  #Return intervals.
  return(interv)
}

#' MPP Segment.
#'
#' Segments an MPP into individual ROIs.
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y, and the final column corresponding to mark.
#' @param m Number of x intervals.
#' @param n Number of y intervals. Will use m if left NULL.
#' @param xlims Vector of minimum and maximum x values to use.
#' @param ylims Vector of minimum and maximum y values to use.
#' @return Segmented ROIs.
#' @export
mppsegment <- function(points, m, n = NULL, xlims = NULL, ylims = NULL){
  #Calculate intervals.
  interv <- intervals(points, m, n, xlims, ylims)
  #Create new point cloud for each segment and add to list.
  return(lapply(1:(m*n), seg, interval = interv, coords = points))
}

#' Create and save plot function used in Save Segments.
#'
#' Create and save plot function used in Save Segments (dev use only).
#'
#' @param points A marked point pattern. Matrix or data frame with first two columns corresponding to x and y, and the final column corresponding to mark.
#' @param loc Location to save.
#' @param titl Plot title.
#' @return Saves segmentation.
#' @export
saveseg <- function(points, loc, titl){
  #Create plot.
  points <- as.data.frame(points)[,1:3]
  originalname <- colnames(points)[3]
  colnames(points) <- c("x", "y", "label")
  
  #Create plot.
  ggplot2::ggplot(data = points, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = label)) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                   axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.75, "cm"),
                   legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16))
  #Save plot.
  ggplot2::ggsave(paste(loc, titl, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 14.342, units = "cm")
}

#' Save Segments.
#'
#' Saves each segmented ROI output from mppsegment().
#'
#' @param segmentation Output of mppsegment().
#' @param location Save location (file path).
#' @return Plot of each ROI.
#' @export
savesegments <- function(segmentation, location){
  #For each found segment, create and save a plot in a specified location.
  for(i in 1:length(segmentation)){
    #Save segment.
    saveseg(points = segmentation[[i]], loc = location, titl = paste("Plot ", i, sep = ""))
  }
}

#' Simulate Photon Count.
#'
#' Generates photon counts from given GP values.
#'
#' @param points A marked point pattern. Matrix or data frame with first three columns corresponding to x, y and GP.
#' @param photons Vector of photon counts for each point in points.
#' @return A marked point pattern. Matrix with six rows: x, y, GP, green photon count, red photon count and GP error.
#' @export
photonsimulate <- function(points, photons){
  #Crop points.
  points <- points[,1:3]
  #Get original column names.
  cols <- colnames(points)
  #Add suggested photon counts.
  points <- cbind(points, photons)
  #Calculate other photon count.
  points <- cbind(points, photons * (1 - points[,3]) / (1 + points[,3]))
  #Replace infinite values.
  points[points[,3] == -1, 5] <- points[points[,3] == -1, 4]
  points[points[,3] == -1, 4] <- 0
  #Get GP errors.
  points <- cbind(points, (1 / sqrt(points[,4] + points[,5])) * (1 + points[,3]))
  #Update column names.
  colnames(points) <- append(cols, c("photons_1", "photons_2", "gp_error"))
  #Return points.
  return(points)
}

#' Plot Marked Point Pattern.
#'
#' Creates a simple plot of a continuous marked point pattern.
#'
#' @param points A marked point pattern. Matrix or data frame with first three columns corresponding to x, y and mark.
#' @param labels Labels of point pattern columns. Will take from points if left NULL.
#' @param range Vector of minimum and maximum marked values to define colour bar.
#' @param includelegend Boolean dictating whether to show legend (colour bar).
#' @return A plot of a marked point pattern.
#' @export
plotmpp <- function(points, labels = NULL, range = NULL, includelegend = TRUE){
  #Convert to data frame.
  points <- as.data.frame(points)
  
  #Check if labels are already included.
  if(is.null(labels)){
    if(ncol(points) > 2){
      labels <- points[,ncol(points)]
    } else {
      stop("No labels given.")
    }
  }
  #Get range if not given.
  if(is.null(range)){
    range <- c(min(labels), max(labels))
  }
  #Choose legend position.
  if(includelegend){
    pos <- "right"
  } else {
    pos <- "none"
  }
  #Update points.
  points <- cbind(points[,1:2], labels)
  colnames(points) <- c("points_x", "points_y", "point_labels")
  #Plot MPP.
  print(
    ggplot2::ggplot(data = points, ggplot2::aes(x = points_x, y = points_y)) +
      ggplot2::geom_point(ggplot2::aes(color = point_labels)) +
      ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01)) +
      ggplot2::scale_color_viridis_c(limits = range) +
      ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
                     axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(), legend.key.height = ggplot2::unit(2.3, "cm"),
                     legend.key.width = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size = 16),
                     legend.position = pos)
  )
}

#' Plot Marked Point Pattern Discrete.
#'
#' Creates a simple plot of a discrete marked point pattern.
#'
#' @param coords A point pattern. Matrix or data frame with two columns corresponding to x and y.
#' @param clusterIndices Cluster indices of each point.
#' @return A plot of the marked point pattern.
#' @export
plotmppdiscrete <- function(coords, clusterIndices){
  numDimensions <- dim(coords)[2]
  if (numDimensions != 2) {
    stop("Data should be 2D")
  }
  numDetections <- dim(coords)[1]
  if (numDetections != length(clusterIndices)) {
    stop("coords and clusterIndices should have the same length")
  }
  uniqueClusters <- unique(clusterIndices)
  numClusters <- sum(uniqueClusters > 0)
  detectionList <- data.frame(x = coords[, 1], y = coords[, 
                                                          2], clusterIndex = clusterIndices)
  if (numClusters > 0) {
    colorMap <- as.character(randomcoloR::distinctColorPalette(numClusters), 
                             luminosity = "dark")
    if (0 %in% uniqueClusters) {
      clusterColors <- c("black", sample(colorMap, numClusters))
    }
    else {
      clusterColors <- sample(colorMap, numClusters)
    }
    ggplot2::ggplot(detectionList, ggplot2::aes(x = x, y = y, color = factor(clusterIndex))) + 
      ggplot2::geom_point() + ggplot2::theme_bw() + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) + 
      ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
      ggplot2::scale_y_continuous(expand = c(0.01,0.01)) +
      ggplot2::scale_color_manual(values = clusterColors) + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                                           axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                                                                           axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                                                                           axis.ticks.y = ggplot2::element_blank(), legend.position = "none")
  }
  else {
    ggplot2::ggplot(detectionList, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point( 
      color = "black") + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                                                              panel.grid.minor = ggplot2::element_blank()) + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                   axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                                                                                                                   axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                                                                                                                   axis.ticks.y = ggplot2::element_blank(), legend.position = "none") +
      ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
      ggplot2::scale_y_continuous(expand = c(0.01,0.01))
  }
}

#All PLASMA Functions.

#' Point Label Analysis of Super-resolved Membrane Attributes
#'
#' Performs P-Check, JOSEPH and MOM. All data should be presented as a data frame with four or five columns: x coordinates, y coordinates, photon count first channel, photon count second channel and (optional, will be calculated if not given) GP. Can input either a single point cloud or a list of point clouds. Each point cloud should have undergone suitable pre-processing (i.e. only include regions) within the cell.
#'
#' @param points A point cloud or list of point clouds (each defined as above).
#' @param location Output location (file path) for all results to be saved.
#' @param radius Cluster radius to use for P-Check and JOSEPH. Will be estimated if left NULL.
#' @param deviance Deviance to use for JOSEPH. Will be estimated if left NULL.
#' @param split Split points to discretise marked point pattern. Will be estimated from Gaussian fitting if left NULL.
#' @param minpoints Minimum number of points per cluster.
#' @param names Vector of names associated with each point pattern. Will use numeric if left NULL.
#' @return All outputs of P-Check, JOSEPH and MOM, saved at given location.
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
      print(plotmppdiscrete(jclust[,1:2], jclust[,4]))
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
      means <- c()
      if(length(unique(jclust[,4])) > 1){
        for(i in unique(jclust[,4])){
          if(i > 0){
            joe <- jclust[jclust[,4] == i,1:2]
            josephvers <- rbind(josephvers, cbind(joe[chull(joe),], i))
            means <- append(means, mean(jclust[jclust[,4] == i, 3]))
          }
        }
        
        #Get viridis colours based on mean GPs.
        mingp <- min(jclust[,3])
        maxgp <- max(jclust[,3])
        scaledmeans <- floor((means - mingp) / (maxgp - mingp) * 1000)
        vircols <- viridis::viridis(1000)
        cols <- vircols[scaledmeans]
        
        #Create plot.
        colnames(josephvers) <- c("x", "y", "v")
        print(ggplot2::ggplot(josephvers, ggplot2::aes(x = x, y = y, fill = factor(v))) +
                ggplot2::geom_polygon() +
                ggplot2::scale_x_continuous(limits = c(minx, maxx), expand = c(0.01, 0.01)) +
                ggplot2::scale_y_continuous(limits = c(miny, maxy), expand = c(0.01, 0.01)) +
                ggplot2::scale_fill_manual(values = cols) +
                ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
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

#' Bulk Analyse function used in PLASMA.
#'
#' Runs each algorithm used in PLASMA (dev use only).
#'
#' @param cloud Point cloud.
#' @param wholeradius Radius used in P-Check and JOSEPH.
#' @param wholedeviance Deviance for JOSEPH.
#' @param split Threshold for discretising MPP.
#' @param minpoints Minimum points per cluster.
#' @return Outputs of each algorithm.
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
    mpphistogramresults <- mom(cloud, location = location, plotname = paste("GP Histogram ", cloudnumber, sep = ""), useggplot = TRUE, outputdata = TRUE)
  } else {
    #Update data for JOSEPH and MPP Histogram results.
    josephresults <- paste("No evidence of separation found (p = ", presults[[1]][2], ")", sep = "")
    mpphistogramresults <- paste("No evidence of separation found (p = ", presults[[1]][2], ")", sep = "")
  }
  
  #Return results.
  return(list("P-Check Results" = presults, "Joseph Results" = josephresults, "MOM Results" = mpphistogramresults))
}

#' Extract Results function used in PLASMA.
#'
#' Extract Results function used in PLASMA (dev use only).
#'
#' @param res Results from bulkanalyse.
#' @param resnum Results number to extract.
#' @return Algorithm result requested.
#' @export
extractresults <- function(res, resnum){
  return(res[[resnum]])
}

#' Plot Histogram Results function used in PLASMA.
#'
#' Plot Histogram Results function used in PLASMA (dev use only).
#'
#' @param i Index.
#' @param bindata Output from bindatum().
#' @param GP_histbins GP values of bins.
#' @param location Save location (file path).
#' @param points Point cloud.
#' @return MPP histogram as per MOM specification.
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

#P-Check.

#' Get Probability Discrete function used in P-Check.
#'
#' Get Probability Discrete function used in P-Check (dev use only).
#'
#' @param points Marked point pattern.
#' @param inradius Logical matrix specifying if points are in radius of each other.
#' @param radtotal Total number of points in radius for each point.
#' @return Probability of non-randomness.
#' @export
getprobability <- function(points, inradius, radtotal){
  #Get probabilities of picking categories within given radius.
  prob <- sum(inradius * (1 - Rfast::Dist(as.double(points[,3]), method = "manhattan")))/radtotal
  return(prob)
}

#' Probability Check function.
#'
#' Performs P-Check, which determines the p-value of non-random partitioning of points by spatial and discrete marked values (probability of domains appearing by random chance).
#'
#' @param points A marked point pattern. Matrix or data frame with first three columns corresponding to x, y and discretised mark.
#' @param number_of_trials Number of permutations to run.
#' @param radius Suspected domain radius (will be estimated if left NULL).
#' @param outputall Output all results (dev use only).
#' @return P-value.
#' @export
pcheck <- function(points, number_of_trials = 100, radius = NULL, outputall = FALSE){
  #Check if points has too many columns.
  points <- points[,c(1,2,ncol(points))]
  #Calculate radius if not given
  if(is.null(radius)){
    radius <- informradius(points)
  }
  #Get inradius.
  inradius <- Rfast::Dist(points[,1:2]) < radius
  #Rename categories.
  categories <- unique(points[,ncol(points)])
  cats <- length(categories)
  number_of_points <- nrow(points)
  if(is.character(categories[1])){
    for(i in 1:cats){
      points[grepl(categories[i], points[,ncol(points)]),ncol(points)] <- i
    }
    categories <- 1:cats
  }
  #Get total points in radius.
  radtotal <- sum(inradius)
  #Calculate original probability.
  originalprobability <- getprobability(points, inradius, radtotal)
  #Randomise the labels and calculate entropy again.
  probs <- sort(do.call(c, lapply(1:number_of_trials, function(i){
    #Randomise labels.
    points[,3] <- sample(points[,3])
    return(getprobability(points, inradius, radtotal))
  })), decreasing = TRUE)
  #Determine the original entropy's percentile.
  for(i in 1:length(probs)){
    if(originalprobability > probs[i]){
      break
    }
  }
  percentile <- i/number_of_trials
  #Output data.
  if(outputall){
    #Return p-values and ARIs as list.
    return(list("Original Weighted Probability and p-value" = c(originalprobability, percentile), "All Probability Values" = probs))
  } else {
    #Return p-values only.
    return(percentile)
  }
}