##---------------------------------------------------------------------------------------------------
# dimemsion-reduce method
#' Reduce dimension of high dimension FCS data
#'
#' We apply dimension reduction to those matrix which has more than paramenters
#' \code{t-SNE},\code{largeVis} and \code{SOM} was used
#'
#' @param object cytosee object
#' @param n_core how many cores you want to use
#' @param sgd_batches input parameter for largeVis (refer largeVis for details)
#' @param tsne_pca whether use pca in tsne
#' @import largeVis
#' @import Rtsne
#' @importFrom destiny DiffusionMap
#' @export
pro_reduce_dim <- function(data,n_core = NULL, sgd_batches = 0.5, tsne_pca=TRUE,method="ALL" ){
    red_dim <-list()
  if(method=="ALL"){
    # run pca
    message("Run pca...")
    for(i in 1:length(data)){
      data[,i]==data[,i]+i*10^-22
    }
    PCA=princomp(data,cor=TRUE)
    red_dim[["PCA"]]<-PCA

    # run largeVis to reduce the dimensions
    message("Run largeVis...")
    suppressMessages(vis <-largeVis(scale(t(data)),sgd_batches=0.8))
    red_dim[["largeVis"]]<-vis

    # run t-SNE to reduce the dimension
    message("Run t-SNE...")
    tsne2d <- Rtsne::Rtsne(data,dims=2, pca=tsne_pca, initial_dims=ncol(data))
    red_dim[["tsne2d"]]<-tsne2d

    # run diffusion map
    message("Run diffusion map")
    diffusion_2d <- DiffusionMap(data)
    red_dim[["diffusion2d"]]<-diffusion_2d
    return(red_dim)
  }
  else if(method=="t-SNE"){

    # run t-SNE to reduce the dimension
    message("Run t-SNE...")
    tsne2d <- Rtsne::Rtsne(data,dims=2, pca=tsne_pca, initial_dims=ncol(data))
    red_dim[["tsne2d"]]<-tsne2d
    return(red_dim)
  }
  else if(method=="PCA"){
    # run pca
    message("Run pca...")
    PCA=princomp(data,cor=TRUE)
    red_dim[["PCA"]]<-PCA
    return(red_dim)
  }
  else if(method=="DiffusionMap"){
    # run diffusion map
    message("Run diffusion map")
    diffusion_2d <- DiffusionMap(data)
    red_dim[["diffusion2d"]]<-diffusion_2d
    return(red_dim)
  }
  else if(method=="LargeVis"){
    # run largeVis to reduce the dimensions
    message("Run largeVis...")
    suppressMessages(vis <-largeVis(scale(t(data)),sgd_batches=0.8))
    red_dim[["largeVis"]]<-vis
    return(red_dim)
  }
}



