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
#' @importFrom  largeVis largeVis
#' @importFrom  Rtsne Rtsne
#' @export
#'
pro_reduce_dim <- function(data,label,sample_n=1000,n_core = NULL, sgd_batches = 0.5, tsne_pca=TRUE ){
  red_dim <-list()
  # run pca

  message("Run pca...")
  PCA=princomp(data,cor=TRUE)
  red_dim[["PCA"]]<-PCA

  smalldata<-data[sample_for_vis(label,sample_n),]

  # run largeVis to reduce the dimensions
  message("Run largeVis...")
  suppressMessages(vis <-largeVis(scale(t(smalldata)),sgd_batches=0.8))
  red_dim[["largeVis"]]<-vis

  # run t-SNE to reduce the dimension
  message("Run t-SNE...")
  tsne2d <- Rtsne::Rtsne(smalldata,dims=2, pca=tsne_pca, initial_dims=ncol(smalldata))
  red_dim[["tsne2d"]]<-tsne2d

  red_dim["red_events"]<-as.data.frame(sample_for_vis(label,sample_n))
  return(red_dim)
}


# sample for dimension reduce
sample_for_vis<-function(label, sample_n = 1000){
  total_len = length(label)
  if(total_len<=sample_n){
    return(1:total_len)
  }
  else{
    label_level = levels(as.factor(label))
    av_num = sample_n / length(label_level)
    return_list=vector()
    for(i in label_level){
      ll = which(label==i)
      if(length(ll)>500){
        return_list = c(return_list,sample(ll,500))
      }
      else if(length(ll)>100){
        return_list = c(return_list,sample(ll,100))
      }
      else {
        return_list = c(return_list,ll)
      }
    }
    return(base::sort(return_list))
  }
}



