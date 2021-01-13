rm(list=ls()); ls()
### libraries #####
require(spatstat); require(spatgraphs); require(sgof); require(dbscan)
set.seed(123)
##############  IMPORT DATA  ####

dd <- read.table("data.csv", sep=";",header=TRUE)

##############  EMST CLUSTERING ########
# Distribution simulation and separation
##############  sub-function  ##########

sim_csr_mst_avg <- function(nsim, lambda, window)  { #function to compute MST edge length under CSR
  set.seed(123)
  avg.dist <- NULL          
  for (i in 1:nsim) {
    sim <- rpoispp(lambda = lambda, win = window)
    mst<-spatgraph(sim ,"MST")        
    if (mst$N>0) { edge.dist=edgeLengths(mst, sim)  
   avg.dist <- c(avg.dist, sum(edge.dist$d)/edge.dist$n)}  
  } 
  return(avg.dist)
}
### End sub-function 

out_of_control <- 3 # result from control chart tool
win <- disc(1,c(0,0))

data_cluster <- dd[dd$lag<=3,c("X","Y")]
data_cluster_pp <- as.ppp(data_cluster[,c("X","Y")], W=win)

#Build MST
mst <- spatgraph(data_cluster_pp ,"MST") 
edge_length <- edgeLengths(mst, data_cluster_pp)
edge_length_sort <- sort(edge_length$d, decreasing = T)

#Simulate distribution
npoint <- data_cluster_pp$n
lambda <- npoint/spatstat::area(win)
start <- Sys.time()
distrSim <- sim_csr_mst_avg(500,lambda, win)
end <- Sys.time()
print(paste("execution time: ",end - start))

#Compute clusters
points <- data_cluster
count <- 0
conditionStop <- round(npoint*.15)
output <- NULL
here <- NULL
pvalue <- NULL

repeat {
  count <- count+1
  pvalue <- c(pvalue, sum(distrSim>edge_length_sort[count])/length(distrSim))
  pvalueAdj <- ifelse(identical(BH(pvalue, .05)$Adjusted.pvalues, numeric(0)), pvalue, BH(pvalue, .05)$Adjusted.pvalues[count])
  rej <- ifelse(pvalueAdj<0.05, 1, 0)
  dbOUT <-  dbscan(cbind(points$X, points$Y), eps=edge_length_sort[count], minPts = 3) # used to label clusters
  numCluster <- as.data.frame(table(dbOUT$cluster)); names(numCluster) <- c("cluster", "Freq")
  numCluster$id <- paste(numCluster$cluster, "db", count)
  here <- NULL
  here <- data.frame(points$X, points$Y, dbOUT$cluster); names(here) <- c("X","Y","cluster")
  here <- merge(here, numCluster, by="cluster")
  here$status <- ifelse(here$Freq<conditionStop, "exit", "repeat")
  output[[count]] <- here[here$status=="exit",]
  points <- here[here$status=="repeat",c("X","Y")]
  if (all(here$status=="exit") | rej==0 | count==length(edge_length_sort)) {output[[count]] <- here; break}
}

clusters <- Reduce(rbind, output) # save the clusters for alpha-shape

#
clustersNN <- clusters[clusters$cluster!=0,]
clustersNN$id <- as.factor(clustersNN$id)
levels(clustersNN$id) <- seq(1,14)

plot(win, main="")
plot(mst, data_cluster_pp, add=T, lty=2,lwd=2,col="lightgray")
for (i in clustersNN$id) {
  plot(spatgraph(clustersNN[clustersNN$id==i,c("X", "Y")] ,"MST"), 
       clustersNN[clustersNN$id==i,c("X", "Y")], add=T, lwd=2 )
}
points(clustersNN$X, clustersNN$Y, pch=as.numeric(clustersNN$id), col=as.numeric(clustersNN$id))
text(0.1878967, 0.3000137 , "A",cex=1.5)
text(0.5557725, 0.364556 , "B",cex=1.5)
text(0.7761832,-0.3669567 , "C",cex=1.5)
text(-0.3679905,-0.4704845, "D",cex=1.5)
