##------------------------------------------------------------------------------------------------------
##cytosee
##cytosee project

#' The cytosee Class
#'
#' The cytosee class is the key of each step of cytosee analysis pipeline. It will records all
#' perfermance and other information. including this dataset, raw data, plot object, analysis, etc.
#' All that is needed to construct a cytosee dataset is a FCS file input.
#'
#' Slots for cytosee
#'
#' @section Slots:
#'  \describe{
#'   \item{\code{raw.data}:}{\code{"matrix"}, The raw exprs data from FCS file }
#'   \item{\code{raw.desc}:}{\code{"list"}, The raw description data from FCS file }
#'   \item{\code{ident}:}{\code{"factor"}, The identified cell types }
#'   \item{\code{projectname}:}{\code{"character"}, The name of your project}
#'   \item{\code{paras_in}:}{\code{"list "}, The input parameters }
#'   \item{\code{dim_red}:}{\code{"list "}, Reduced dimemsion data }
#'   \item{\code{pro.pre}:}{\code{"list "}, The pre-processing information }
#'   \item{\code{pro.clust}:}{\code{"list "}, clustering data info }
#'   \item{\code{outermethod}:}{\code{"list "}, clustering data info }
#'  }
#'
#' @name cytosee
#' @rdname cytosee
#' @aliases cytosee-class
#' @exportClass cytosee


cytosee <- methods::setClass("cytosee",slots =
                      c(autoLabel="logical",
                        channel.use = "numeric",
                        clust="list",
                        ClusterID="list",
                        clust_method="character",
                        CL_label="list",
                        dim.red="list",
                        dmt="list",
                        event.use = "numeric",
                        fcs.data = "ANY",
                        filename = "character",
                        ID2CL="list",
                        label="ANY",
                        mst="ANY",
                        outermethod="list",
                        preprocess="character",
                        projectname="character",
                        transform_method = "character",
                        transform_paras = "character",
                        version = "ANY")
)



# Hide cytosee data
setMethod("show","cytosee",function(object){
  cat("An object of ",class(object)," class for the analysis of cytometry data")
  invisible(NULL)
})
