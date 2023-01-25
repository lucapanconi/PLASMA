#P-Check.
#Get Probability Discrete function used in Entropy Check.
getprobability <- function(points, inradius, radtotal){
  #Get probabilities of picking categories within given radius.
  prob <- sum(inradius * (1 - Rfast::Dist(as.double(points[,3]), method = "manhattan")))/radtotal
  return(prob)
}

#Probability Check function.
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