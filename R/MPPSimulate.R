#Marked Point Pattern Simulator.
#Toroidal Distance Matrix function used in MPP Simulate.
toroidaldist <- function(points, domain_max_x, domain_max_y){
  #Compute x and y differences.
  x_d <- do.call(cbind, lapply(1:nrow(points), function(i){
    x_d <- abs(points[i,1] - points[,1])
    x_d[x_d > domain_max_x] <-  domain_max_x * 2 - x_d[x_d > domain_max_x]
    return(x_d)
  }))
  y_d <- do.call(cbind, lapply(1:nrow(points), function(i){
    y_d <- abs(points[i,2] - points[,2])
    y_d[y_d > domain_max_y] <-  domain_max_y * 2 - y_d[y_d > domain_max_y]
    return(y_d)
  }))
  #Create toroidal coordinates matrix.
  return(sqrt(x_d ** 2 + y_d ** 2))
}

#Calculate Getis L Values function used in MPP Simulate.
calculate_L_values <- function(distance_matrix, number_of_points, domain_size_x, domain_size_y, search_radius){
  return(do.call(c, lapply(1:number_of_points, function(i){
    return(sqrt(((domain_size_x * domain_size_y) / (pi * number_of_points)) * (sum(rep(1,number_of_points)[distance_matrix[i,] <= search_radius])-1)))
  })))
}

#Velocity function used in MPP Simulate.
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

#MPP Simulate.
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
  cluster_index <- RSMLM::clusterTomato(points, tomato_radius, tomato_threshold)
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

#Segment function used in MPP Segment.
seg <- function(x, interval, coords){
  #Take the given interval and extract only points in that interval.
  return(coords[coords[,1] >= interval[x,1] & coords[,1] < interval[x,2] &
                coords[,2] >= interval[x,3] & coords[,2] < interval[x,4],])
}

#Intervals function used in MPP Segment.
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

#MPP Segment.
mppsegment <- function(points, m, n = NULL, xlims = NULL, ylims = NULL){
  #Calculate intervals.
  interv <- intervals(points, m, n, xlims, ylims)
  #Create new point cloud for each segment and add to list.
  return(lapply(1:(m*n), seg, interval = interv, coords = points))
}

#Create and save plot function used in Save Segments.
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

#Save Segments.
savesegments <- function(segmentation, location){
  #For each found segment, create and save a plot in a specified location.
  for(i in 1:length(segmentation)){
    #Save segment.
    saveseg(points = segmentation[[i]], loc = location, titl = paste("Plot ", i, sep = ""))
  }
}

#Photon Count.
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

#Plot Marked Point Pattern.
plotmpp <- function(points, labels = NULL, range = NULL, includelegend = TRUE){
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

#Plot Marked Point Pattern Discrete.
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
                                                                                                            panel.grid.minor = ggplot2::element_blank()) + theme(axis.title.x = ggplot2::element_blank(), 
                                                                                                                                                        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                                                                                                                                                        axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                                                                                                                                                        axis.ticks.y = ggplot2::element_blank(), legend.position = "none") +
      ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
      ggplot2::scale_y_continuous(expand = c(0.01,0.01))
  }
}