#Justification of Separation by Employed Persistent Homology
#' JOSEPH
#' 
#' Outputs a partition of the data set, clustering by distance and marks.
#' 
#' @param points Numeric point pattern.
#' @param deviance Maximum difference in marked values between roots of clusters and any other point.
#' @param radius Search radius defining neighbourhood of each point.
#' @return Partition of data set.
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
    
    #Extract neighbours with lower filtration values that i.
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

#' Filter JOSEPH
#' 
#' Filters raw JOSEPH results either by minimum number of points per cluster, range of GP values per cluster or whether dispersion is displayed at the given radius value in Ripley’s H function.
#' 
#' @param josephresults Results from JOSEPH.
#' @param mode Mode of filtering. Must be a vector containing any combination of the characters "dispersion", "minpoints" and "range".
#' @param minpoints Minimum number of points allowed in cluster.
#' @param range Range of marked values allowed in cluster.
#' @return Partition of data set.
#' @export
filterjoseph <- function(josephresults, mode = "dispersion", minpoints = NULL, range = NULL, radius = NULL){
  #Determine number of clusters.
  number_of_clusters <- max(unique(josephresults[,4]))
  #Check each possible mode.
  if("minpoints" %in% mode){
    #Check minpoints given.
    if(is.null(minpoints)){
      stop("No minimum given. Add argument minpoints = n, where n is the desired minimum number of points per cluster.")
    } else {
      #Iterate over all clusters.
      for(i in 1:number_of_clusters){
        #Find points corresponding to that cluster.
        cluster_indicator <- josephresults[,4] == i
        if(sum(cluster_indicator) < minpoints){
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
      #Iterate over all clusters.
      for(i in 1:number_of_clusters){
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
    #Iterate over all clusters.
    for(i in 1:number_of_clusters){
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

#' Split JOSEPH
#' 
#' Splits results from JOSEPH into separate data sets based on the splits found by inform persistence. Takes the average label of each cluster and sorts it into one of several categories, then displays each category separately.
#' 
#' @param josephresults Results from JOSEPH.
#' @param splits Split points for discretisation, will estimate by Gaussian Mixture Modelling if omitted.
#' @return Partition of data set.
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
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size=1),
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

#' JOSEPH Cluster Summary
#' 
#' Calculates summary statistics from JOSEPH results.
#' 
#' @param josephresults Results from JOSEPH.
#' @return Summary statistics of JOSEPH results.
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