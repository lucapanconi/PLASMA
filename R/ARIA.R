#Adjusted Rand Index Analysis.
#Discretise.
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
                       legend.title = ggplot2::element_text(size = 20), legend.text = ggplot2::element_text(size = 16))
    )
  }
  #Return new points.
  return(points)
}

#Inform Cluster Persistence.
informclusterpersistence <- function(points, radius = NULL){
  #Calculate number of neighbours within radius.
  if(is.null(radius)){
    radius <- informradius(points)
  }
  #Calculate number of neighbours within radius for each point.
  numberneighbours <- as.numeric(rowSums(Rfast::Dist(points[,1:2]) <= radius)) - 1
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

#Get ARI function used in ARIA.
getari <- function(points, clusterpersistence = NULL, clusterradius = NULL, createplot = FALSE, location = "", plotname = "",
                   clusterwhole = NULL){
  #Check if points has too many columns.
  points <- points[,c(1,2,ncol(points))]
  #Calculate cluster radius if not given
  if(is.null(clusterradius)){
    clusterradius <- informradius(points)
  }
  #Calculate cluster persistence if not given.
  if(is.null(clusterpersistence)){
    clusterpersistence <- informclusterpersistence(points, radius = clusterradius)
  }
  #Convert point cloud to data frame if not already
  if(is.matrix(points)){
    points <- as.data.frame(points)
  }
  #Determine number of categories
  categories <- sort(unique(points[,3]))
  number_of_categories <- length(categories)
  #ToMATo needs at least two categories to run
  if(number_of_categories == 1){
    #Create standard ARI matrix.
    ARIs <- matrix(0L, ncol = 2, nrow = 1)
    #Update values (no partition, so ARI will be 1 regardless)
    ARIs[1,1] <- "No Partition"
    ARIs[1,2] <- 1
  } else{
    #Rename discrete categories to numbers
    category_numbers <- 1:number_of_categories
    
    #Calculate number of points
    number_of_points <- nrow(points)
    
    #Performs ToMATo on the entire data set unless whole clustering is already specified
    if(is.null(clusterwhole)){
      clusterwhole <- RSMLM::clusterTomato(points[,1:2], clusterradius, clusterpersistence) 
    }
    
    #Devise each of the possible partitions
    partitions <- partitions::listParts(number_of_categories)
    number_of_partitions <- length(partitions)
    ARIs <- matrix(0L, ncol = 2, nrow = number_of_partitions - 1)
    
    #For each partition (except the first as we require at least two distinct sets), pass through each set in the partition,
    #extract only the points belonging to the categories in that set and cluster them, combine their clusterings and compare
    #to the whole using the ARI
    #For each partition, extract the partition
    for(i in 2:number_of_partitions){
      partition <- partitions[[i]]
      cluster_index_individual <- c()
      point_index <- c()
      
      #Set empty partition name
      partition_name <- ""
      
      #Store all points and their subsets for plotting
      all_points <- data.frame()
      sub_set_number <- 0
      
      #For each extracted partition, extract each subset
      for(j in 1:length(partition)){
        sub_set <- partition[[j]]
        sub_set_number <- sub_set_number + 1
        #For each subset, extract the points in its categories
        sub_set_points <- matrix(0L, ncol = 2)
        #Update partition name
        partition_name <- paste(partition_name, "{", sep = "")
        for(k in sub_set){
          #Determine subset points
          cat <- categories[k]
          sub_set_points <- rbind(sub_set_points, as.matrix(points[points[,3] == cat, 1:2]))
          point_index <- append(point_index, (1:number_of_points)[points[,3] == cat])
          #Add category to partition name
          partition_name <- paste(partition_name, cat, sep = " ")
        }
        #Update partition name
        partition_name <- paste(partition_name, " } ", sep = "")
        #Remove empty row
        sub_set_points <- as.matrix(sub_set_points[-1,])
        #Store points and their subset number
        all_points <- rbind(all_points, cbind(sub_set_points, rep(sub_set_number, nrow(sub_set_points))))
        #Record number of clusters found so far
        total_number_of_clusters <- sum(unique(cluster_index_individual) >= 0)
        #If there is at least two points, perform ToMATo on that category
        if(nrow(sub_set_points) == 1){
          cluster_index_individual <- append(cluster_index_individual, 0)
        } else{
          #Perform ToMATo on that category
          cluster_index_individual <- append(cluster_index_individual, RSMLM::clusterTomato(sub_set_points, clusterradius, clusterpersistence) + total_number_of_clusters)
        }
      }
      #Update column names
      colnames(all_points) <- c("points_x","points_y","point_labels")
      #Plot and save all points
      if(createplot){
        #Create plot
        print(
          p <- ggplot2::ggplot() +
            ggplot2::geom_point(data = all_points, ggplot2::aes(x = points_x, y = points_y, colour = as.factor(point_labels)), size = 2) +
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
        #Save plot
        ggplot2::ggsave(paste(location, plotname, " ", partition_name, ".pdf", sep = ""), device = "pdf", width = 17.342, height = 17.342, units = "cm")
      }
      
      #Combine individual cluster index with point index and reorder
      cluster_index_individuals <- as.matrix(as.data.frame(cbind(cluster_index_individual, point_index))[order(point_index),])[,1]
      
      #Calculate adjusted rand index
      ARIs[i - 1, 1] <- partition_name
      ARIs[i - 1, 2] <- mclust::adjustedRandIndex(clusterwhole, cluster_index_individuals)
    }
  }
  #Update column names
  colnames(ARIs) <- c("Partition", "ARI")
  rownames(ARIs) <- 1:nrow(ARIs)
  
  #Return all ARIs
  return(ARIs)
}

#ARIA.
aria <- function(points, number_of_trials = 100, clusterpersistence = NULL, clusterradius = NULL, outputall = FALSE){
  #Check if points has too many columns.
  points <- points[,c(1,2,ncol(points))]
  #Calculate cluster radius if not given
  if(is.null(clusterradius)){
    clusterradius <- informradius(points)
  }
  #Calculate cluster persistence if not given.
  if(is.null(clusterpersistence)){
    clusterpersistence <- informclusterpersistence(points, radius = clusterradius)
  }
  if(is.null(clusterpersistence)){
    cat("Parameter estimation failed, not enough data to determine radius and/or persistence. Returning NULL.\n")
    return(NULL)
  }
  #Calculate whole clustering.
  clusterwhole <- RSMLM::clusterTomato(points[,1:2], clusterradius, clusterpersistence) 
  #Determine original ARIs.
  original_ARIs <- getari(points, clusterpersistence, clusterradius, clusterwhole = clusterwhole)
  #Shuffle the original point labels.
  original_labels <- points[,3]
  ARIs <- original_ARIs[,1]
  #Randomise the labels and perform the ARI test again.
  for(i in 1:number_of_trials){
    #Randomise labels.
    points[,3] <- sample(points[,3], replace = FALSE)
    current_ARIs <- getari(points, clusterpersistence, clusterradius, clusterwhole = clusterwhole)
    ARIs <- cbind(ARIs, current_ARIs[,2])
  }
  #Sort each row and determine the original ARIs percentile.
  percentiles <- original_ARIs
  percentage <- 1 / number_of_trials
  sorted_ARIs <- data.frame()
  for(i in 1:nrow(ARIs)){
    original_value <- original_ARIs[i,2]
    values <- sort(ARIs[i, 2:(number_of_trials + 1)], decreasing = TRUE)
    sorted_ARIs <- rbind(sorted_ARIs, values)
    for(j in 1:number_of_trials){
      if(values[j] <= original_value){
        break
      }
    }
    percentiles[i,2] <- percentage * j
  }
  #Add p-values.
  percentiles <- cbind(original_ARIs, percentiles[,2])
  #Update row and column names for percentiles.
  colnames(percentiles) <- c("Partition", "ARI", "Estimated p-value")
  rownames(percentiles) <- 1:nrow(percentiles)
  #Output data.
  if(outputall){
    #Update column names for all ARI values.
    colnames(sorted_ARIs) <- 1:number_of_trials
    #Return p-values and ARIs as list.
    return(list("ARI and p-value" = percentiles, "All ARI Values" = sorted_ARIs))
  } else {
    #Return p-values only.
    return(percentiles)
  }
}