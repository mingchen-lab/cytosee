## Calculate adaboost model
##----------------------------------------------------------------------------------
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



# densitycut flowSOM flowmeans rphenograph
##----------------------------------------------------------------------------------
#' @title The flowMeans algorithm
#' @description Finds a good fit to the data using k-means clustering algorithm.
#' Then merges the adjacent dense spherical clusters to find non-spherical clusters.
#'
#' @name cytosee_flowMeans
#' @aliases cytosee_flowMeans
#'
#' @param object An object of \code{cytosee} class.
#' @param varNames A character vector specifying the variables (columns) to be included in clustering.
#'  When it is left unspecified, all the variables will be used.
#' @param MaxN Maximum number of clusters. If set to NA (default) the value will be estimated automatically.
#' @param NumC Number of clusters. If set to NA (default) the value will be estimated automatically.
#' @param iter.max The maximum number of iterations allowed.
#' @param nstart The number of random sets used for initialization.
#' @param Mahalanobis Boolean value. If TRUE (default) mahalanobis distance will be used. Otherwised, euclidean distance will be used.
#' @param Standardize Boolean value. If TRUE (default) the data will be transformed to the [0,1] interval.
#' @param Update String value. If set to "Mahalanobis" the distance function will be updated at each merging iteration with recalculating mahalanobis distances.
#'  If set to "Mean" the distance matrix will be updated after each merging step with averaging. If set to "None" the distance matrix will not be updated.
#' @param MaxCovN Maximum number of points, used for calculating the covariance. If set to NA (default), all the points will be used.
#' @param MaxKernN Maximum number of points, used for counting the modes using kernel density estimation. If set to NA (default), all the points will be used.
#' @param addNoise Boolean value. Determines if uniform noise must be added to the data to prevent singularity issues or not.
#' @param OrthagonalResiduals Boolean value, indicates if the residuals must be transformed to orthagonal distance or not.
#'
#'
#' @importFrom flowMeans flowMeans
#' @export
#'
#'

cytosee_flowMeans <- function(object, varNames=NULL, MaxN = NA, NumC = NA, iter.max = 50, nstart = 10,
                              Mahalanobis = TRUE, Standardize = TRUE, Update = "Mahalanobis", OrthagonalResiduals=TRUE,
                              MaxCovN=NA, MaxKernN=NA, addNoise=TRUE){
  data=object@fcs.data[object@event.use,][,object@channel.use]
  clust <- flowMeans::flowMeans(data,MaxN=MaxN,NumC=NumC,varNames=varNames,iter.max = iter.max, nstart = nstart,
                                Mahalanobis = Mahalanobis, Standardize = Standardize, Update = Update, OrthagonalResiduals=OrthagonalResiduals,
                                MaxCovN=MaxCovN, MaxKernN=MaxKernN, addNoise=addNoise)
  clust <- clust@Label
  object@ClusterID <- as.data.frame(clust)
  return(object)
}



##----------------------------------------------------------------------------------
#' @title RphenoGraph clustering
#'
#' @description R implementation of the PhenoGraph algorithm
#' A simple R implementation of the [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm,
#' which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing
#' phenotypic similarities between cells by calclating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities
#' using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph.
#'
#' @param object an object of \code{cytosee} class.
#' @param k integer; number of nearest neighbours (default:30).
#'
#' @return a list contains an igraph graph object for \code{graph_from_data_frame} and a communities object, the operations of this class contains:
#' \item{print}{returns the communities object itself, invisibly.}
#' \item{length}{returns an integer scalar.}
#' \item{sizes}{returns a numeric vector.}
#' \item{membership}{returns a numeric vector, one number for each vertex in the graph that was the input of the community detection.}
#' \item{modularity}{returns a numeric scalar.}
#' \item{algorithm}{returns a character scalar.}
#' \item{crossing}{returns a logical vector.}
#' \item{is_hierarchical}{returns a logical scalar.}
#' \item{merges}{returns a two-column numeric matrix.}
#' \item{cut_at}{returns a numeric vector, the membership vector of the vertices.}
#' \item{as.dendrogram}{returns a dendrogram object.}
#' \item{show_trace}{returns a character vector.}
#' \item{code_len}{returns a numeric scalar for communities found with the InfoMAP method and NULL for other methods.}
#' \item{plot}{for communities objects returns NULL, invisibly.}
#'
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. Cell, 2015.
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' Rphenograph_out <- Rphenograph(data, k = 45)
#' modularity(Rphenograph_out[[2]])
#' membership(Rphenograph_out[[2]])
#' iris_unique$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))
#' ggplot(iris_unique, aes(x=Sepal.Length, y=Sepal.Width, col=Species, shape=phenograph_cluster)) + geom_point(size = 3)+theme_bw()
#' @import igraph
#' @export
cytosee_phenograph <- function(object, K = 30){
  data <- object@fcs.data[object@event.use,][,object@channel.use]
  clust <- Rphenograph(data,k = K)
  clust <- igraph::membership(clust[[2]])
  object@ClusterID <- as.data.frame(as.matrix(clust))
  return(object)
}


##----------------------------------------------------------------------------------
#' @title The densityCut algorithm
#'
#' @description densityCut first roughly estimates the densities
#'of data points from a K-nearest neighbour graph and then refines the densities via a random
#'walk. A cluster consists of points falling into the basin of attraction of an estimated mode of the
#'underlining density function. A post-processing step merges clusters and generates a hierarchical
#'cluster tree. The number of clusters is selected from the most stable clustering in the hierarchical
#'cluster tree. Experimental results on ten synthetic benchmark datasets and two microarray gene
#'expression datasets demonstrate that densityCut performs better than state-of-the-art algorithms
#'for clustering biological datasets.
#'
#' @name cytosee_DensityCut
#' @aliases cytosee_DensityCut
#' @export
#'
#' @param object an object of \code{cytosee} class.
#' @param X A data matrix (columns are features and rows are data points)
#' @param K A integer to specify the number of neighbours in building the Knn graph.
#' Default to \eqn{K=\log_2(N)}, where N is the number of data points
#' @param knn.index An N*K data matrix for the nearest neighbour indices
#' @param knn.dist An N*K data matrix for the nearest neighbour distances
#'
#' @param threshold A number between 0 and 1 specifying the saliency index to cut the tree.
#' If not specified, it is selecting by stability analysis of the clustering tree
#'
#' @param V The initial density vector of length N
#' @param D The dimensionality of data
#' @param G A sparse Knn graph, reseaved for extension
#'
#' @param    The damping factor between 0 and 1, default to 0.90
#' @param nu The saliency index in merging trees, default to \eqn{seq(0.0, 1.0, by=0.05)}
#' @param adjust Lotical, whether to ajdust valley height or not
#' @param alpha The damping factor between 0 and 1, default to 0.90
#' @param maxit The maximum number of iteration allowed in density refinement, default to 50
#' @param eps The threshold in density refinement, default to 1e-5
#'
#' @param col A vector of clours
#' @param show.plot Logical, whether to draw clustering results
#' @param show.tip.label Logical, whether to draw the tip labels of trees
#'
#' @param debug Logical, whether to print debug information
#' @param xlab Logical, whether to show the xlab
#' @param text subplot label
#' @param ... Reserved for extension
cytosee_DensityCut <- function(object, K, knn.index, knn.dist, V, D, G, threshold,
                               alpha=0.90,adjust=TRUE, maxit=50, eps=1e-5,
                               col = NULL,debug=FALSE, xlab=TRUE, text=NULL, ...){
  X<- object@fcs.data[object@event.use,][,object@channel.use]
  clust <- DensityCut(X, K=K, knn.index, knn.dist, V, D, G, threshold=threshold,
                      alpha=alpha, adjust=TRUE, maxit=maxit, eps=eps,
                      col = NULL, show.plot=FALSE,show.tip.label=FALSE,
                      debug=debug, xlab=xlab, text=text, ...)$cluster
  object@ClusterID <- as.data.frame(clust)
  return(object)
}




##----------------------------------------------------------------------------------
#' @title Run the flowSOM algorithm
#'
#' @description
#' Method to run general FlowSOM workflow. Will scale the data and uses consensus meta-clustering.
#' @param object an object of \code{cytosee} class.
#' @param K Maximum number of clusters to try out.
#' @param nClus Exact number of clusters to use. If not NULL, max will be ignored.
#' @param method Clustering method to use, given as a string. Options are metaClustering_consensus,metaClustering_hclust, metaClustering_kmeans,metaClustering_som
#'
#' @name cytosee_flowSOM
#' @aliases cytosee_flowSOM
#' @import FlowSOM
#' @export
cytosee_flowSOM=function(object,K,nClus=NULL,method="metaClustering_consensus"){
  File <- flowFrame(as.matrix(object@fcs.data[object@event.use,]))
  ff <- File
  colsToUse <- cyto@channel.use
  fSOM <- ReadInput(ff,compensate = FALSE,transform = FALSE, scale = TRUE)
  fSOM <- BuildSOM(fSOM,colsToUse = colsToUse,silent = TRUE)
  metaClustering <- MetaClustering(data=fSOM$map$codes,max=K,method=method,nClus=nClus)
  result <- matrix()
  for(i in 1:length(fSOM$map$mapping[,1])){
    result[i] <- metaClustering[fSOM$map$mapping[i,1]]
  }
  object@ClusterID <- as.data.frame(result)
  return(object)
}


##----------------------------------------------------------------------------------
#' SamSPECTRAL image based funciton for fcs data clustering.
#'
#' SamSPECTRAL first builds the communities to sample the data points. Then, it builds a graph and after weighting
#' the edges of the graph by conductance computation, it is passed to a classic spectral clustering algorithm to find
#' the spectral clusters. The last stage of SamSPECTRAL is to combine the spectral clusters.
#' The resulting "connected components" estimate biological cell populations in the data sample.
#'
#' @param object an object of \code{cytosee} class.
#' @param sigmal A scaling parameter that determines the "resolution" in the spectral clustering stage.
#'  By increasing it, more spectral clusters are identified. This can be useful when "small" population are aimed.
#'  See the user manual for a suggestion on how to set this parameter using the eigenvalue curve.
#' @param separation.factor This threshold controls to what extend clusters should be combined or kept separate.Normally,
#'  an appropriate value will fall in range 0.3-2.
#' @param number.of.clusters The default value is "NA" which leads to computing the number of spectral clusters automatically,
#'  otherwise it can be a vector of integers each of which determines the number of spectral clusters.
#'  The output will contain a clustering resulting from each value.
#'
#' @return A object of cytosee with clustering result.
#'
#' @name cytosee_SamSPECTRAL
#' @aliases cytosee_SamSPECTRAL
#' @import SamSPECTRAL
#' @export
#'
cytosee_SamSPECTRAL<-function(object,sigma,separation.factor,number.of.clusters="NA",precision=6,stabilizer=100,k.for_kmeans=NA){
  data <- object@fcs.data[object@event.use,]
  data <- as.matrix(data)
  result <- SamSPECTRAL::SamSPECTRAL(
                                  data.points=data,
                                  dimensions = object@channel.use,
                                  normal.sigma = sigma,
                                  separation.factor=separation.factor,
                                  number.of.clusters=number.of.clusters,
                                  precision=precision,
                                  stabilizer = stabilizer,
                                  k.for_kmeans=k.for_kmeans

  )
  object@ClusterID <- as.data.frame(result)
  return(object)
}
