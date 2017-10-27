#############################################
# Data Generation
#############################################
set.seed(11)

# pick 2 main colors to act as cluster centers
col_label <- c("green","yellow")
r <- c(0, 240)
g <- c(158, 228)
b <- c(115, 66)
col_palette <- cbind.data.frame(r,g,b,col_label)
size <- c(10,20)

# simulate similar color values around the 2 main colors
d <- NULL
n <- 500
for (i in 1:length(col_label))  {
  sd <- size[i]
  r <- rnorm(n=n, mean=col_palette$r[i], sd=sd)
  g <- rnorm(n=n, mean=col_palette$g[i], sd=sd)
  b <- rnorm(n=n, mean=col_palette$b[i], sd=sd)
  label <- col_palette$col_label[i]

  # consolidate into dataset
  rgb <- cbind.data.frame(r,g,b,label)
  d <- rbind.data.frame(d,rgb)
}

# remove out-of-range RGB values
rgb_out <- which(d[,1:3]<0 | d[,1:3]>255, arr.ind = T)
rgb_out <- unique(rgb_out[,1])
d <- d[-rgb_out,]

# add hex values, more convenient to feed into R plotting commands
d$hex <- rgb(d$r, d$g, d$b, max=255)

# plot and save chosen colors
chosen_dim <- c(2,3)
plot(d[,chosen_dim],
     ylab = "Blue Value",
     xlab = "Green Value",
     col=d$hex, pch=16,
     xlim=c(130,255), ylim=c(0,150))

#############################################
# Data Analysis
#############################################

# install and load package
require(SOMbrero)

# set number of iterations
niter <- 500

# train SOM
grid_dim <- 8
set.seed(11)
my.som <- trainSOM(d[,1:3], dimension=c(grid_dim, grid_dim),
                   nb.save=niter/2, maxit=niter,
                   radius.type="letremy")

# plot convergence
plot(my.som, what="energy")

# plot map on each iteration
npt <- grid_dim^2
steps <- my.som$backup$steps

for (s in 1:length(steps))  {

  # get neuron coordinates
  iterate_map <- my.som$backup$prototypes[[s]][,chosen_dim]

  # plot
  plot(d[,chosen_dim], col=d$hex,
       pch=16, cex=.5,
       main=paste('Iteration ', steps[s]),
       ylim=c(0,150),
       xlim=c(130,255))
  points(iterate_map)

  for (pt in 1:npt) {

    # identify grid neighbors
    pt_neighbor <- c(pt-1,pt+1,
                     pt+grid_dim,pt-grid_dim)

    # remove out-of-range neighbors
    rm_out <- which(pt_neighbor < 1 | pt_neighbor > npt)

    # remove criss-cross links
    if ((pt %% grid_dim) == 0) {
      rm_crisscross <- which(pt_neighbor == (pt+1))
    }
    if ((pt %% grid_dim) == 1) {
      rm_crisscross <- which(pt_neighbor == (pt-1))
    }
    pt_neighbors <- pt_neighbor[-c(rm_out,rm_crisscross)]

    for(i in pt_neighbors)  {
      lines(x <- c(iterate_map[pt,1],iterate_map[i,1]),
            y <- c(iterate_map[pt,2],iterate_map[i,2]))
    }
  }
}


# plot counts in each neuron
plot(my.som, what="obs", type="hitmap")
table(my.som$clustering)

# radar plot for each neuron
plot(my.som, what="obs", type="radar",
     key.loc=c(-0.5,5), mar=c(0,10,2,0))

# plot distance between neurons
plot(my.som, what="prototypes", type="smooth.dist")

# label neurons with freq of color occurrence
# make dummy matrix
dummy1 <- as.numeric(d$label==chosen_col[1])
dummy2 <- as.numeric(d$label==chosen_col[2])
dummy_labels <- cbind(dummy1, dummy2)
colnames(dummy_labels) <- chosen_col
# plot
plot(my.som, what="add",
     type="words", variable=dummy_labels)

# label neurons based on color prediction
plot(my.som, what="add",
     type="color", variable=dummy_labels[,1])


####################################################
# Clustering of Neurons
####################################################

my.clusters <- superClass(my.som, k=2)
summary(my.clusters)

plot(my.clusters, plot.var=FALSE)
plot(my.clusters, type='grid', plot.legend=TRUE)
plot(my.clusters, type="hitmap", plot.legend=TRUE)


####################################################
# Making Predictions
####################################################

library(kohonen)

# train model
my.som2 <- xyf(as.matrix(d[,1:3]), as.factor(d$label),
               grid=somgrid(grid_dim, grid_dim, "rectangular"))

# generate predictions
som.predictions <- predict(my.som2)

# confusion matrix
table(som.predictions$predictions[[2]], d$label)

# distance between data point and best matching unit (BMU)
# possible proxy for uncertainty
boxplot(my.som2$distances ~ d$label)

# prediction for each neuron
som.predictions$unit.predictions

