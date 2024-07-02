# Point Label Analysis of Super-resolved Mark Attributes Membrane
Point Label Analysis of Super-resolved Mark Attributes, or PLASMA, is a data analysis package written in R for simulating, interpreting and visualising marked point patterns, with an emphasis on those derived from single molecule localisation microscopy (SMLM). The algorithms included can be run entirely through [**RStudio**](https://www.rstudio.com/products/rstudio/).

## Installation
The package can be installed using pre-existing functions in the [**devtools**](https://cran.r-project.org/web/packages/devtools/index.html) package. It is highly recommended that you use **RStudio**. If you do not have **devtools** installed, enter the following line into the command prompt:
```{r}
install.packages("devtools")
```
Then, simply enter the line below into the command prompt:
```{r}
devtools::install_github("lucapanconi/PLASMA")
```

## Simulation
### Simulating Domains
In the example code below, we simulate a marked point pattern of 500 points with a cluster radius of 50 units. Each point is assigned a random value from one of two Normal distributions, depending on whether it falls in a domain (category) or the background.
```{r}
#Implement PLASMA library.
library(PLASMA)

#Include mclust (you only need to run this line once).
require(mclust)

#Simulate point pattern.
points <- mppsimulate(number_of_points = 500, search_radius = 50, number_of_categories = 1, category_mean = c(1), category_sd = c(0.1), background_mean = -1, background_sd = 0.3)

#Plot marked point pattern.
plotmpp(points[,1:3])
```

<img src="https://github.com/lucapanconi/PLASMA/blob/master/Simulation1.png" width="801">

By default, the probability of a point being assigned to a domain is 0.7. For a full breakdown of all parameters in the simulator, run the line *?mmpsimulate*.

Note that the [**mclust**](https://cran.r-project.org/web/packages/mclust/index.html) package is required for many of the functions in PLASMA and can sometimes be temperamental when installing. The easiest fix for this is to include the line *require(mclust)* at the beginning of your code.

### Multiple categories.
We can increase the number of categories - that is, distributions of values in domains - by tuning the parameters. For example, in the simulation below, we generate three different category types and assign points to each of them with probability 0.3. These domains have means 1, 2 and 5 and standard deviations 0.1, 0.3 and 1, respectively. Each simulation will display at least one domain of each type of category when possible.
```{r}
#Simulate point pattern with multiple categories.
points <- mppsimulate(number_of_points = 500, search_radius = 50, number_of_categories = 3, category_mean = c(1, 2, 5), category_sd = c(0.1, 0.3, 1), probability_of_picking_categories = c(0.3, 0.3, 0.3), background_mean = -1, background_sd = 0.3)

#Plot marked point pattern.
plotmpp(points[,1:3])
```

<img src="https://github.com/lucapanconi/PLASMA/blob/master/Simulation2.png" width="801">

## Usage
The main features of PLASMA are the algorithms used for determining whether there is some degree of co-occurrence of points with similar values. There are three main algorithms, known as P-Check, JOSEPH and MPP Histogram. All of these functions can be accessed individually or used in batch by the main algorithm, *plasma*, as seen below.

```{r}
#Perform all algorithms in PLASMA.
plasma(points, location = "C:/Users/user_name/Desktop/output_folder/")
```

Here, *points* can be given as a single marked point pattern (usually represented as a data frame) or a list of point patterns. Each point pattern must contain at least four columns, representing the x coordinates, y coordinates, photon counts in first channel (lower wavelength) and photon counts in second channel, respectively. All results will be output to the given location. For more details on each constituent function, see below.

### Probability Check (Checking for Existence of Domains)
The first tool introduced in the PLASMA repertoire is known as Probability Check, or P-Check. First and foremost, this algorithm aims to determine if a discrete point pattern expresses non-random heterogeneity or, more succinctly, if domains are present within the data. The method itself operates on the premise that data sets exhibiting non-random colocalisation of identically-labelled points will express, on average, neighbourhoods which are predominantly homogeneous. We initialise the algorithm with a search radius from which we can determine the neighbourhood of each point. The search radius will be estimated from the [**Ripley's H Function**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2726315/) if not provided.

```{r}
#P-Check requires a discrete point pattern.
discrete_points <- discretise(points)

#Perform P-Check on point cloud.
pcheck(discrete_points, radius = 50)
```

This outputs an estimate of the *p* value representing the probability that co-occurrence of identically-labelled points has occurred by random chance. The smaller the value, the more likely it is that domains are in your data set.

### Justification of Separation by Employed Persistent Homology (Clustering Points in Domains)
Once the existence of domains had been established by P-Check, we can utilise advanced cluster analysis techniques to identify their locations and coverage. To do this, we introduce another novel technique, known as Justification of Separation by Employed Persistent Homology (JOSEPH). Building on existing topological data analysis (TDA) methods, this algorithm probes for collections of points which are both spatially correlated and express marks lying within a quantifiable range of their mean value. Essentially, JOSEPH acts as a clustering algorithm for marked point patterns.

```{r}
#Perform JOSEPH on a point cloud.
josephresults <- joseph(points[,1:3], radius = 50, deviance = 2)

#Plot cluster results.
plotmppdiscrete(josephresults[,1:2], josephresults[,4])
```

<img src="https://github.com/lucapanconi/PLASMA/blob/master/JosephResults.png" width="801">

Unlike PLASMA, P-Check and JOSEPH only require point patterns with three columns, representing the *(x,y)* coordinates and the mark for each point. This makes them applicable to any generic marked point pattern data. JOSEPH does not differentiate outliers by default. To do this, use the *filterjoseph* function as seen below.

```{r}
#Filter results.
filteredresults <- filterjoseph(josephresults)

#Plot filtered results.
plotmppdiscrete(filteredresults[,1:2], filteredresults[,4])
```

<img src="https://github.com/lucapanconi/PLASMA/blob/master/FilteredResults.png" width="801">

For a full breakdown of all filtering methods, run the line *?filterjoseph*.

To acquire summary statistics of each cluster, use the function *josephclustersummary*.

```{r}
#Get summary statistics.
summarystats <- josephclustersummary(filteredresults)
```

This includes the number of points, area, density, average mark, mark standard deviation and difference from global mean of each cluster.

### MPP Histogram (Visualising Generalised Polarisation Values)

In the context of SMLM, we may wish to probe the degree of membrane order using solvatochromatic probes (polarity-sensitive dyes). The raw data from this acquisition will take the form of a point cloud, with each point assigned a photon count from two channels, therefore yielding two marks. We can combine these two measurements into one metric, known as the Generalised Polarisation (GP) value, which represents the degree of membrane order at each localisation. PLASMA contains functions for both calculating and displaying GP values.

```{r}
#Get GP Values.
gpvalues <- gpcalculate(points)

#Plot distribution of GP Values as histogram, with density plot overlaid.
gphistogram(gpvalues, density = TRUE)
```

<img src="https://github.com/lucapanconi/PLASMA/blob/master/GPHistogram.png" width="801">

Here, the *points* data frame has four columns, representing the x coordinates, y coordinates, photon counts in first channel (lower wavelength) and photon counts in second channel, respectively. Visualising GP Values can be difficult in dense point patterns. As such, we can bin photon counts into a 2D histogram and use this to calculate a GP image at a higher resolution.

```{r}
#Create MPP Histogram.
mom(points, location = "C:/Users/user_name/Desktop/output_folder/")
```

This will output a PDF, containing the histograms of the photons counts in each channel and the GP values, to the given location. For more information on Generalised Polarisation, see [**Owen 2010**](https://pubmed.ncbi.nlm.nih.gov/19937746/).

## License
Licensed under GNU General Public License v3.0.