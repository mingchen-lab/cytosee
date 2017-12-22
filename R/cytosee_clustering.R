#' Opens \code{cytosee} results in an interactive session in a web browser.
#' @param object an object of \code{cytosee} class
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize
#' all precomputed clusterings.
#'
#' @name cytosee_clustering
#' @aliases cytosee_clustering





##### flowSOM #####

#' @import FlowSOM
#' @export


cytosee_flowSOM=function(File,K,colsToUse){
  ff <- File
  fSOM <- ReadInput(ff,compensate = FALSE,transform = FALSE, scale = TRUE)
  fSOM <- BuildSOM(fSOM,colsToUse = colsToUse)
  metaClustering <- metaClustering_consensus(data=fSOM$map$codes,k=K)
  result=matrix()
  for(i in 1:length(fSOM$map$mapping[,1])){
    result[i]=metaClustering[fSOM$map$mapping[i,1]]
  }
  return(result)
}

##---------------------------------------------------------------------------------------------------
# Function for calculate the center of clusters
# pro_center_calc for building the mst
pro_center_calc <- function(data,label) {
  centers <- c()
  is.na(label) <-which(label == 0)
  for(i in c(1:max(label, na.rm=TRUE))) {
    obs <- which(label == i)
    if (length(obs) > 1) {
      centers <- rbind(centers,colMeans(data[obs,,drop=FALSE]))
      label[obs] <- nrow(centers)
    } else {
      is.na(label) <- obs
    }

  }
  return(centers)
}

estkTW <- function(dataset) {

  p <- ncol(dataset)
  n <- nrow(dataset)

  # compute Tracy-Widom bound
  x <- scale(dataset)
  muTW <- (sqrt(n - 1) + sqrt(p))^2
  sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
  #  sigmaHatNaive <- tmult(x)  # x left-multiplied by its transpose
  sigmaHatNaive <- x %*% t(x)  # x left-multiplied by its transpose
  bd <- 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution

  # compute eigenvalues and return the amount which falls above the bound
  evals <- eigen(sigmaHatNaive, symmetric = TRUE, only.values = TRUE)$values
  k <- 0
  for (i in 1:length(evals)) {
    if (evals[i] > bd) {
      k <- k + 1
    }
  }
  return(k)
}

pic_mad <- function(data){
  m_mad<-apply(data,1,mad)
  data <- data[sort(m_mad,index.return=T,decreasing = T)$ix[1:5000],]
  return(data)
}

##---------------------------------------------------------------------------------------------------
# Function for calculate the distance matrix
# pro_center_calc for building the mst
pro_calc_distance <- function(data) {
  distMat <- list()
  distMat[[1]]<- as.matrix(dist(data,method = "euclidean"))
  distMat[[2]]<- as.matrix(dist(data,method = "maximum"))
  distMat[[3]]<- as.matrix(dist(data,method = "manhattan"))
  distMat[[4]]<- as.matrix(dist(data,method = "minkowski"))
  return(distMat)
}



##---------------------------------------------------------------------------------------------------
# calculate_elbow, this script comes from FlowSOM

#' Generate the best label by consensus clustering
#'
#' We use consensus clustering method to make sure the result is reliable
#'
#' @param object cytosee object
#' @param MaxC The max number to esitimate the k
#' @param n_core how many cores you want to use
#' @export

PRO_consensus <- function(object , MaxC = MaxC , n_cores = n_cores){

  if(method ==    "metaClustering_consensus"){
    results <- consensus(data,max,...)
    res <- rep(0,max)
    res[1] <- SSE(data,rep(1,nrow(data)))
    for(i in 2:max){
      c <- results[[i]]$consensusClass
      res[i] <- SSE(data,c)
    }
  } else {
    method <- get(method)
    res <- rep(0,max)
    for(i in 1:max){
      c <- method(data, k=i,...)
      res[i] <- SSE(data,c)
    }
  }

  for(i in 2:(max-1)){
    res[i] <- (1-smooth)*res[i]+(smooth/2)*res[i-1]+(smooth/2)*res[i+1]
  }

  if(plot) plot(1:max, res, type="b", xlab="Number of Clusters",
                ylab="Within groups sum of squares")
  findElbow(res)
}

findElbow <- function(data){
  n <- length(data)
  data <- as.data.frame(cbind(1:n,data))
  colnames(data) <- c("X","Y")

  min_r <- Inf
  optimal <- 1
  for(i in 2:(n-1)){
    f1 <- stats::lm(Y~X,data[1:(i-1),])
    f2 <- stats::lm(Y~X,data[i:n,])
    r <- sum(abs(c(f1$residuals,f2$residuals)))
    if(r < min_r){
      min_r <- r
      optimal <-i
    }
  }
  return(optimal)
}

SSE <- function(data,clustering){
  if(class(clustering)!= "numeric")
    clustering <- as.numeric(as.factor(clustering))
  c_wss <- 0
  for(j in seq_along(clustering)){
    if(sum(clustering==j) > 1){
      c_wss <- c_wss + (nrow(data[clustering==j,,drop=FALSE])-1)*
        sum(apply(data[clustering==j,,drop=FALSE],2,stats::var))
    }
  }
  return(c_wss)
}


################### temp method use flowmeans to calculate the clustering info

clusteringMethod <- function(dat, MaxN=20){
  result <- flowMeans(dat, MaxN=MaxN )
  result <- list(result@Label, result@Labels, result@Mats)
  return(result)
}



########### consensus with consensusPlus
## require ada and caret

# sampling cells randomly
downsamleRD<- function(data,K_cell = 4000) {
  sample_index <- sample(1:nrow(data),K_cell)
}

oneVsAll <- function(X,Y,FUN,...) {
  models <- lapply(unique(Y), function(x) {
    name <- as.character(x)
    .Target <- factor(ifelse(Y==name,name,'other'), levels=c(name, 'other'))
    dat <- data.frame(.Target, X)
    model <- FUN(.Target~., data=dat, ...)
    return(model)
  })
  names(models) <- unique(Y)
  info <- list(X=X, Y=Y, classes=unique(Y))
  out <- list(models=models, info=info)
  class(out) <- 'oneVsAll'
  return(out)
}

predict.oneVsAll <- function(object, newX=object$info$X, ...) {
  stopifnot(class(object)=='oneVsAll')
  lapply(object$models, function(x) {
    predict(x, newX, ...)
  })
}

classify <- function(dat) {
  out <- dat/rowSums(dat)
  out$Class <- apply(dat, 1, function(x) names(dat)[which.max(x)])
  out
}


#' culculate the adjusted Rand Index for the return class with the gd class
#' @param x the return class from clustering method
#' @param y the standard label
adjustedRandIndex <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}





# Rclustpp densitycut flowSOM flowmeans rphenograph

#' @export
cytosee_Rclusterpp <- function(object, K=8, n_cores = 1){
  Rclusterpp::Rclusterpp.setThreads(threads = n_cores)
  clust <- Rclusterpp::Rclusterpp.hclust(object)
  cluster_id <- cutree(clust,K)
  return(cluster_id)
}


#' @importFrom flowMeans flowMeans
#' @export
cytosee_flowMeans <- function(object, NumC=8, MaxN = MaxN ){
  clust <- flowMeans::flowMeans(object,MaxN=MaxN,NumC=NumC)
  clust <- clust@Label
  return(clust)
}

#' @import Rphenograph
#' @import igraph
#' @export
cytosee_phenograph <- function(object, K = 30){
  clust <- cytofkit::Rphenograph(object,k = K)
  clust <- igraph::membership(clust)
  return(clust)
}

#' @export
cytosee_densitycut <- function(object, NumC= 30){
  clust <- DensityCut(object,show.plot = F,K=NumC)$cluster
  return(clust)
}

##### flowSOM #####

#' @import FlowSOM
#' @export
cytosee_flowSOM=function(File,K,colsToUse){
  ff <- File
  fSOM <- ReadInput(ff,compensate = FALSE,transform = FALSE, scale = TRUE)
  fSOM <- BuildSOM(fSOM,colsToUse = colsToUse)
  metaClustering <- metaClustering_consensus(data=fSOM$map$codes,k=K)
  result=matrix()
  for(i in 1:length(fSOM$map$mapping[,1])){
    result[i]=metaClustering[fSOM$map$mapping[i,1]]
  }
  return(result)
}
