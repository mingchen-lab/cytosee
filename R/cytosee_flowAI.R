# Detects anomalies in a time series using Cyclic hybrid ESD (C-H-ESD).
#
# Anomaly Detection Using Cyclic Hybrid ESD Test ----GM
#
# A technique for detecting anomalies in univariate time
# series where the input is a series of observations.
# @name AnomalyDetection
# @param x Time series as a column data frame, list, or vector,
#  where the column consists of the observations.
# @param max_anoms Maximum number of anomalies that C-H-ESD will
# detect as a percentage of the data. The value can be from 0 to 1.
# @param direction Directionality of the anomalies to be detected.
# Options are: \code{'pos' | 'neg' | 'both'}.
# @param alpha The level of statistical significance with which
# to accept or reject anomalies.
# @param use_decomp If set to \code{'FALSE'} it gives the possibility
# to detect outliers with the generalized ESD method on the orginal data.
# By default is set to \code{'TRUE'} and time series decomposition
# is performed before the analysis.
# @param period Defines the number of observations in a single
# period, and used during seasonal decomposition.
# @param e_value Add an additional column to the anoms output
# containing the expected value.
# @param verbose Additionally printing for debugging.

# @details
# @return The returned value is a list with the following components.
# @return \item{anoms}{Data frame containing index, decomposition components, and
# optionally expected values.}
# @return One can save \code{anoms} to a file in the following fashion:
# \code{write.csv(<return list name>[["anoms"]], file=<filename>)}
# @references Rosner, B., (May 1983), "Percentage Points for a
# Generalized ESD Many-Outlier Procedure", Technometrics, 25(2),
# pp. 165-172.
# @export


anomaly_detection = function(x, max_anoms=0.49, direction='both', alpha=0.01, use_decomp = TRUE, period=1, verbose = FALSE){

  idNOzero <- which(x != 0)
  x <- x[idNOzero]

  # Check for supported inputs types
  if(is.vector(x) && is.numeric(x)) {
    x <- ts(x, frequency = period)
  } else if(is.ts(x)) {
  } else {
    stop("data must be a time series object or a vector that holds numeric values.")
  }

  # Handle NAs
  if (length(rle(is.na(c(NA,x,NA)))$values)>3){
    stop("Data contains non-leading NAs. We suggest replacing NAs with interpolated values (see na.approx in Zoo package).")
  } else {
    x <- na.omit(x)
  }

  # Sanity check all input parameterss
  if(max_anoms > .49){
    stop(paste("max_anoms must be less than 50% of the data points (max_anoms =", round(max_anoms*length(x), 0), " data_points =", length(x),")."))
  }
  if(!direction %in% c('pos', 'neg', 'both')){
    stop("direction options are: pos | neg | both.")
  }
  if(!(0.01 <= alpha || alpha <= 0.1)){
    print("Warning: alpha is the statistical significance level, and is usually between 0.01 and 0.1")
  }
  if(is.null(period)){
    stop("Period must be set to the number of data points in a single period")
  }

  ############## -- Main analysis: Perform C-H-ESD -- #################
  # -- Step 1: Decompose data. This will return two more components: trend and cycle
  if(use_decomp){
    x_cf <- cffilter(x)
    #med_t <- trunc(median(x_cf$trend))
    med_t <- trunc(median(x))
    sign_n <- sign(x_cf$trend - med_t)
    sign_n[which(sign_n == 0)] <-1
    # add the absolute values of the cycle component to the absolute values of the centered trend component. The signs are then added again
    x_2 <- as.vector(trunc(abs(x - med_t) + abs(x_cf$cycle)) * sign_n)
  } else {
    x_2 <- as.vector(x - median(x))
  }

  anomaly_direction = switch(direction,
                             "pos" = data.frame(one_tail=TRUE, upper_tail=TRUE), # upper-tail only (positive going anomalies)
                             "neg" = data.frame(one_tail=TRUE, upper_tail=FALSE), # lower-tail only (negative going anomalies)
                             "both" = data.frame(one_tail=FALSE, upper_tail=TRUE)) # Both tails. Tail direction is not actually used.


  n <- length(x_2)
  data_det <- data.frame(index = idNOzero, values = x_2, or_values = x)
  # Maximum number of outliers that C-H-ESD can detect (e.g. 49% of data)
  max_outliers <- trunc(n*max_anoms)
  func_ma <- match.fun(median)
  func_sigma <- match.fun(mad)
  R_idx <- 1L:max_outliers
  num_anoms <- 0L
  one_tail <- anomaly_direction$one_tail
  upper_tail <- anomaly_direction$upper_tail
  # Compute test statistic until r=max_outliers values have been
  # removed from the sample.
  for (i in 1L:max_outliers){
    if(verbose) message(paste(i,"/", max_outliers,"completed"))

    if(one_tail){
      if(upper_tail){
        ares <- data_det[[2L]] - func_ma(data_det[[2L]])
      } else {
        ares <- func_ma(data_det[[2L]]) - data_det[[2L]]
      }
    } else {
      ares = abs(data_det[[2L]] - func_ma(data_det[[2L]]))
    }

    # protect against constant time series
    data_sigma <- func_sigma(data_det[[3L]])
    # the standard deviation has to be calculated from the orginal
    # distribution because otherwise it would be affected too much
    # by the cycle component
    if(data_sigma == 0)
      break

    ares <- ares/data_sigma
    R <- max(ares)

    temp_max_idx <- which(ares == R)[1L]

    R_idx[i] <- data_det[[1L]][temp_max_idx]

    data_det <- data_det[-which(data_det[[1L]] == R_idx[i]), ]

    ## Compute critical value.
    if(one_tail){
      p <- 1 - alpha/(n-i+1)
    } else {
      p <- 1 - alpha/(2*(n-i+1))
    }

    t <- qt(p,(n-i-1L))
    lam <- t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))

    if(R > lam)
      num_anoms <- i
  }

  if(num_anoms > 0) {
    R_idx <- R_idx[1L:num_anoms]
    all_data <- data.frame(index = idNOzero, anoms = x)
    anoms_data <- subset(all_data, (all_data[[1]] %in% R_idx))
  } else {
    anoms_data <- NULL
  }
  return (list(anoms = anoms_data, num_obs = n))
}

#' Automatic quality control of flow cytometry data.
#'
#' For a set of FCS files, \emph{flow_auto_qc} performs a complete and automatic
#' quality control. It consists in the detection and removal of
#' anomalies by checking three properties of flow cytometry: 1) flow rate,
#' 2) signal acquisition, 3) dynamic range.
#'
#' @param fcsfiles It can be a character vector with the filenames of the FCS files,
#' a flowSet or a flowFrame.
#' @param remove_from Select from which of the three steps the anomalies have to be
#' excluded in the high quality FCS file. The default option \code{"all"} removes the
#' anomalies from all the three steps. Alternatively, you can use:
#' \code{"FR_FS", "FR_FM", "FS_FM", "FR", "FS", "FM"}, to remove the anomalies only
#' on a subset of the steps where \emph{FR} stands for the flow rate, \emph{FS} stands
#' for signal acquisition and \emph{FM} stands for dynamic range.
#' @param output Set it to 1 to return a flowFrame or a flowSet with high quality events only.
#' Set it to 2 to return a flowFrame or a flowSet with an additional
#' parameter where the low quality events have a value higher than 10,000.
#' Set it to 3 to return  a list with the IDs of low quality cells. Set it to any other
#' value if no R object has to be returned. Default is \code{1}.
#' @param timeCh Character string corresponding to the name of the Time Channel
#' in the set of FCS files. By default is \code{NULL} and the name is retrieved
#' automatically.
#' @param second_fractionFR The fraction of a second that is used to split
#' the time channel in order to recreate the flow rate. Set it to
#' \code{"timestep"} if you wish to recreate the flow rate at the maximum
#' resolution allowed by the flow cytometry instrument. Usually, the timestep
#' corresponds to 0.01, however, to shorten the running time of the analysis the
#' fraction used by default is 0.1, corresponding to 1/10 of a second.
#' @param alphaFR The level of statistical significance used to
#' accept anomalies detected by the ESD method. The default value is \code{0.01}.
#' @param decompFR Logical indicating whether the flow rate should be decomposed
#' in the trend and cyclical components. Default is \code{TRUE} and the ESD
#' outlier detection will be executed on the trend component penalized by the
#' magnitude of the cyclical component. If it is \code{FALSE} the ESD outlier
#' detection will be executed on the original flow rate.
#' @param ChRemoveFS Add a character vector with the names or name portions
#' of the channels that you want to exclude from the signal acquisition check.
#' The default option, \code{c("FSC", "SSC")}, excludes the scatter parameters.
#' If you want to include all the parameters in the analysis use \code{NULL}.
#' @param outlierFS logical indicating whether outliers have to be removed
#' before the changepoint detection of the signal acquisition check.
#' The default is \code{FALSE}.
#' @param pen_valueFS The value of the penalty for the changepoint detection
#' algorithm. This can be a numeric value or text giving the formula to use;
#' for instance, you can use the character string \code{"1.5*log(n)"},
#' where n indicates the number of cells in the FCS file. The higher the
#' penalty value the less strict is the detection of the anomalies.
#' The default is \code{200}.
#' @param max_cptFS The maximum number of changepoints that can be detected
#' for each channel. The default is \code{3}.
#' @param ChFM A character vector that indicates which channels need to include
#' for the dynamic range check. The default option is \code{NULL} and
#' with it all the channels are selected for the analysis.
#' @param sideFM Select whether the dynamic range check has to be executed on
#' both limits, the upper limit or the lower limit. Use one of the options:
#' \code{"both", "upper", "lower"}. The default is \code{"both"}.
#' @param neg_valuesFM Scalar indicating the method to use for the removal of the
#' anomalies from the lower limit of the dynamic range. Use \code{1} to remove
#' negative outliers or use \code{2} to truncate the negative values to the cut-off
#' indicated in the FCS file.
#' @param html_report Suffix to be added to the FCS filename to name the HTML
#' report of the quality control. The default is \code{"_QC"}. If you do not
#' want to generate a report use \code{FALSE}.
#' @param mini_report Name for the TXT file containing the percentage of
#' anomalies detected in the set of FCS files
#' analyzed. The default is \code{"_QCmini"}. If you do not want to generate
#' the mini report use \code{FALSE}.
#' @param fcs_QC Suffix to be added for the filename of the new FCS
#' containing a new parameter where the low quality events only have a value
#' higher than 10,000. The default is \code{"_QC"}.
#' If you do not want to generate the high quality FCS file use \code{FALSE}.
#' @param fcs_highQ Suffix to be added for the filename of the new FCS
#' containing only the events that passed the quality control. The default
#' is \code{FALSE} and hence the high quality FCS file is not generated.
#' @param fcs_lowQ Suffix to be added for the filename of the new FCS
#' containing only the events that did not pass the quality control. The default
#' is \code{FALSE} and hence the low quality FCS file is not generated.
#' @param folder_results Character string used to name the directory that
#' contains the results. The default is \code{"resultsQC"}. If you intend
#' to return the results in the working directory use \code{FALSE}.
#' @return A complete quality control is performed on flow cytometry data in FCS
#' format. By default the analysis returns:
#'
#' 1. a flowFrame or flowSet object containing new FCS files with only high quality events
#'
#' and a directory named \var{resultsQC} containing:
#'
#' 1. a set of new FCS files with a new parameter to gate out the low quality events a value larger than 10,000 is assigned to them only,
#'
#' 2. a set of HTML reports, one for each FCS file, that include graphs and table indicating where the anomalies were detected,
#'
#' 3. a single TXT file reporting the percentage of events removed in each FCS file.
#'
#'
#' @import plyr
#' @import knitr
#' @import reshape2
#' @importFrom changepoint cpt.meanvar
#' @importFrom scales pretty_breaks
#' @importFrom graphics hist legend lines
#' @importFrom methods as is new
#' @importFrom stats convolve frequency is.ts mad median na.omit qt runif ts tsp
#' @importFrom utils write.table
#' @export
cytosee_flowAI <- function(fcsfiles, remove_from = "all", output = 1,
                         timeCh = NULL, second_fractionFR = 0.1, alphaFR = 0.01, decompFR = TRUE,
                         ChRemoveFS = c("FSC", "SSC"), outlierFS = FALSE, pen_valueFS = 200,
                         max_cptFS = 3, ChFM = NULL, sideFM = "both", neg_valuesFM = 1,
                         html_report = "_QC", mini_report = "QCmini", fcs_QC = "_QC", fcs_highQ = FALSE,
                         fcs_lowQ = FALSE, folder_results = "resultsQC") {

  ## load the data
  if( is.character(fcsfiles) ){
    set <- read.flowSet(files = fcsfiles)
    names <- fcsfiles
  }else if( class(fcsfiles) == "flowSet"){
    set <- fcsfiles
    names <- flowCore::sampleNames(fcsfiles)
  }else if( class(fcsfiles) == "flowFrame" ){
    set <- as(fcsfiles,"flowSet")
    names <- fcsfiles@description$GUID
  }else{
    stop("As first argument, use a flowSet or a character vector with the path of the FCS files")
  }

  N_cell_set <- flow_set_qc(set)
  area.color <- rep("red", length(set))

  if (missing(timeCh) || is.null(timeCh)) {
    timeCh <- findTimeChannel(set[[1]])
  }
  if (is.null(timeCh)) {
    warning("Impossible to retrieve the time channel automatically. The quality control can only be performed on signal acquisition and dynamic range.", call. =FALSE)
  }

  # in some cases, especially if the FCS file has been modified, there
  # could be more than one slots for the Timestep parameter. the first in
  # numerical order should be the original value.
  word <- which(grepl("TIMESTEP", names(set[[1]]@description),
                      ignore.case = TRUE))
  timestep <- as.numeric(set[[1]]@description[[word[1]]])
  if( !length(timestep) ){
    warning("The TIMESTEP keyword was not found and hence it was set to 0.01. Graphs labels indicating time might not be correct", call. =FALSE)
    timestep <- 0.01
  }
  if( second_fractionFR == "timestep" ){
    second_fractionFR <- timestep
  }else if( second_fractionFR < timestep ){
    stop("The argument second_fractionFR must be greater or equal to timestep.", call. =FALSE)
  }

  if(folder_results != FALSE){
    folder_results <- strip.sep(folder_results)
    dir.create(folder_results, showWarnings = FALSE)
    folder_results <- paste0(folder_results, .Platform$file.sep)
  } else { folder_results <- ""}

  out <- list()

  for (i in 1:length(set)) {
    filename_ext <- basename(description(set[[i]])$FILENAME)
    filename <- sub("^([^.]*).*", "\\1", filename_ext)

    if (html_report != FALSE) {
      reportfile <- paste0(folder_results,filename, html_report, ".html")
    }
    if (mini_report != FALSE) {
      minireport <-  paste0(folder_results, mini_report, ".txt")
      if(!file.exists(minireport)){
        write.table(t(c("Name file", "n. of events", "% anomalies", "analysis from",
                        "% anomalies flow Rate",  "% anomalies Signal",  "% anomalies Margins")),
                    minireport, sep="\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
      }
    }
    if (fcs_QC != FALSE) {
      QC.fcs.file <- paste0(folder_results, filename, fcs_QC, ".fcs")
    }
    if (fcs_highQ != FALSE) {
      good.fcs.file <- paste0(folder_results, filename, fcs_highQ, ".fcs")
    }
    if (fcs_lowQ != FALSE) {
      bad.fcs.file <- paste0(folder_results,filename, fcs_lowQ, ".fcs")
    }
    cat(paste0("Quality control for the file: ", filename, "\n"))
    # select different color for the analyzed FCS in the set plot
    area <- area.color
    area[i] <- "blue"

    # check the time channel of the file
    if (!is.null(timeCh)) {
      if (length(unique(exprs(set[[i]])[, timeCh])) == 1){
        cat("The time channel contains a single value. It cannot be used to recreate the flow rate. \n")
        warning(paste0("The time channel in ", filename_ext, " contains a single value. It cannot be used to recreate the flow rate. \n"), call. =FALSE)
        TimeChCheck <- "single_value"
      }else{
        TimeChCheck <- NULL
      }
    }else{
      TimeChCheck <- "NoTime"
    }

    # get the size of the bins
    FSbinSize <- min(max(1, ceiling(nrow(set[[1]])/100)), 500)
    # order events in the FCS file if a proper Time channel is present
    if (is.null(TimeChCheck)) {
      ordFCS <- ord_fcs_time(set[[i]], timeCh)
    }else{
      ordFCS <- set[[i]]
    }

    origin_cellIDs <- 1:nrow(ordFCS)

    ### Describe here the arguments for the functions of the flow Rate and Flow Signal
    FR_bin_arg <- list( second_fraction = second_fractionFR, timeCh = timeCh,
                        timestep = timestep)
    FR_QC_arg <- list( alpha = alphaFR, use_decomp = decompFR)
    FS_bin_arg <- list( binSize = FSbinSize, timeCh = timeCh, timestep = timestep, TimeChCheck = TimeChCheck)
    FS_QC_arg <- list( ChannelRemove = ChRemoveFS, pen_valueFS, max_cptFS, outlierFS )
    FM_QC_arg <- list( margin_channels = ChFM , side= sideFM, neg_values = neg_valuesFM)

    #### The actual analysis is performed here
    if (is.null(TimeChCheck)) {
      FlowRateData <- do.call(flow_rate_bin, c(ordFCS, FR_bin_arg ))
      FlowRateQC <- do.call(flow_rate_check_a, c(ordFCS, list(FlowRateData), FR_QC_arg))
    }else{
      FlowRateQC<-list()
      FlowRateQC$goodCellIDs <- origin_cellIDs
      FlowRateQC$res_fr_QC$badPerc <- 0
    }
    FlowSignalData <- do.call(flow_signal_bin_a, c(ordFCS,FS_bin_arg))
    FlowSignalQC <- do.call(flow_signal_check_a, c(ordFCS,list(FlowSignalData),FS_QC_arg))
    FlowMarginQC <- do.call(flow_margin_check_a, c(ordFCS, FM_QC_arg))

    ####selection of good cells
    if(remove_from == "all"){
      goodCellIDs <- intersect(FlowRateQC$goodCellIDs, intersect(FlowSignalQC$goodCellIDs, FlowMarginQC$goodCellIDs))
      analysis <- "Flow Rate, Flow Signal and Flow Margin"
    }else if(remove_from == "FR_FS"){
      goodCellIDs <- intersect(FlowRateQC$goodCellIDs, FlowSignalQC$goodCellIDs)
      analysis <- "Flow Rate and Flow Signal"
    }else if(remove_from == "FR_FM"){
      goodCellIDs <- intersect(FlowRateQC$goodCellIDs, FlowMarginQC$goodCellIDs)
      analysis <- "Flow Rate and Flow Margin"
    }else if(remove_from == "FS_FM"){
      goodCellIDs <- intersect(FlowSignalQC$goodCellIDs, FlowMarginQC$goodCellIDs)
      analysis <- "Flow Signal and Flow Margin"
    }else if(remove_from == "FR"){
      goodCellIDs <- FlowRateQC$goodCellIDs
      analysis <- "Flow Rate"
    }else if(remove_from == "FS"){
      goodCellIDs <- FlowSignalQC$goodCellIDs
      analysis <- "Flow Signal"
    }else if(remove_from == "FM"){
      goodCellIDs <- FlowMarginQC$goodCellIDs
      analysis <- "Flow Margin"
    }

    badCellIDs <- setdiff(origin_cellIDs, goodCellIDs)
    totalBadPerc <- round(length(badCellIDs)/length(origin_cellIDs), 4)
    sub_exprs <- exprs(ordFCS)
    params <- parameters(ordFCS)
    keyval <- keyword(ordFCS)
    if (fcs_highQ != FALSE || output == 1) {
      goodfcs <- flowFrame(exprs = sub_exprs[goodCellIDs, ], parameters = params, description = keyval)
      if (fcs_highQ != FALSE) {suppressWarnings(write.FCS(goodfcs, good.fcs.file)) }
    }
    if (fcs_QC != FALSE || output == 2 ){
      QCvector <- FlowSignalData$cellBinID[,"binID"]
      if(length(QCvector) > 9000) QCvector <- runif(length(QCvector), min=1, max=9000)
      QCvector[badCellIDs] <- runif(length(badCellIDs), min=10000, max=20000)
      newFCS <- addQC_a(QCvector, remove_from, sub_exprs, params, keyval)
      if (fcs_QC != FALSE){ suppressWarnings(write.FCS(newFCS, QC.fcs.file)) }
    }
    if (length(badCellIDs) > 0 && fcs_lowQ != FALSE) {
      badfcs <- flowFrame(exprs = sub_exprs[badCellIDs, ],parameters = params,description = keyval)
      suppressWarnings(write.FCS(badfcs, bad.fcs.file))
    }
    if (mini_report != FALSE) {
      write.table(t(c(filename, as.integer(dim(set[[i]])[1]),
                      totalBadPerc * 100, analysis, FlowRateQC$res_fr_QC$badPerc * 100,
                      FlowSignalQC$Perc_bad_cells$badPerc_tot * 100,
                      FlowMarginQC$badPerc * 100)), minireport, sep="\t",
                  append=TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
    }
    if (html_report != FALSE) {
      h_FS_graph <- round(0.4 * (ncol(ordFCS)),1)
      if (!is.null(ChRemoveFS)){
        ChannelRemovedFS <- as.character(grep(paste(ChRemoveFS, collapse="|"),
                                              ordFCS@parameters$name, value = TRUE))
      }
      template_path <- system.file("rmd","autoQC_report.Rmd", package='flowAI')
      knit2html(template_path, output = reportfile, force_v1 = TRUE, quiet = TRUE)
      #   file.remove("autoQC_report.md")
      #   unlink("figure", recursive = TRUE)
    }
    if(output == 1){
      out <- c(out, goodfcs)
    }else if ( output == 2){
      out <- c(out, newFCS)
    }else if( output == 3 ){
      out[[i]] <- badCellIDs
      names(out)[i] <- filename
    }
  }
  if( output == 1 || output == 2){
    if(length(out) == 1){ return( out[[1]] )
    }else{
      OutSet <- as(out, "flowSet")
      flowCore::sampleNames(OutSet) <- names
      return( OutSet ) }
  }
  if( output == 3 ){ return(out) }
}

### Christiano-Fitzgerald filter
cffilter <- function(x,pl=NULL,pu=NULL,root=FALSE,drift=FALSE,
                     type=c("asymmetric","symmetric","fixed","baxter-king",
                            "trigonometric"),nfix=NULL,theta=1)
{

  type = match.arg(type)
  if(is.null(root)) root <- FALSE
  if(is.null(drift)) drift <- FALSE
  if(is.null(theta)) theta <- 1
  if(is.null(type)) type <- "asymmetric"

  if(is.ts(x))
    freq=frequency(x)
  else
    freq=1

  if(is.null(pl))
  {
    if(freq > 1)
      pl=trunc(freq*1.5)
    else
      pl=2
  }

  if(is.null(pu))
    pu=trunc(freq*8)

  if(is.null(nfix))
    nfix = freq*3

  nq=length(theta)-1;
  b=2*pi/pl;
  a=2*pi/pu;

  xname=deparse(substitute(x))
  xo = x
  x = as.matrix(x)
  n = nrow(x)
  nvars = ncol(x)

  if(n < 5)
    warning("# of observations < 5")

  if(n < (2*nq+1))
    stop("# of observations must be at least 2*q+1")

  if(pu <= pl)
    stop("pu must be larger than pl")

  if(pl < 2)
  {
    warning("pl less than 2 , reset to 2")
    pl = 2
  }

  if(root != 0 && root != 1)
    stop("root must be 0 or 1")

  if(drift<0 || drift > 1)
    stop("drift must be 0 or 1")

  if((type == "fixed" || type == "baxter-king") && nfix < 1)
    stop("fixed lag length must be >= 1")

  if(type == "fixed" & nfix < nq)
    stop("fixed lag length must be >= q")

  if((type == "fixed" || type == "baxter-king") && nfix >= n/2)
    stop("fixed lag length must be < n/2")

  if(type == "trigonometric" && (n-2*floor(n/2)) != 0)
    stop("trigonometric regressions only available for even n")

  theta = as.matrix(theta)
  m1 = nrow(theta)
  m2 = ncol(theta)
  if(m1 > m2)
    th=theta
  else
    th=t(theta)

  ##   compute g(theta)
  ##   [g(1) g(2) .... g(2*nq+1)] correspond to [c(q),c(q-1),...,c(1),
  ##                                        c(0),c(1),...,c(q-1),c(q)]
  ##   cc = [c(0),c(1),...,c(q)]
  ## ?? thp=flipud(th)
  g=convolve(th,th,type="open")
  cc = g[(nq+1):(2*nq+1)]
  ##   compute "ideal" Bs
  j = 1:(2*n)
  B = as.matrix(c((b-a)/pi, (sin(j*b)-sin(j*a))/(j*pi)))
  ##    compute R using closed form integral solution
  R = matrix(0,n,1)
  if(nq > 0)
  {
    R0 = B[1]*cc[1] + 2*t(B[2:(nq+1)])*cc[2:(nq+1)]
    R[1] = pi*R0
    for(i in 2:n)
    {
      dj = Bge(i-2,nq,B,cc)
      R[i] = R[i-1] - dj
    }
  }
  else
  {
    R0 = B[1]*cc[1]
    R[1] = pi*R0;
    for(j in 2:n)
    {
      dj = 2*pi*B[j-1]*cc[1];
      R[j] = R[j-1] - dj;
    }
  }

  AA = matrix(0,n,n)

  ###  asymmetric filter

  if(type == "asymmetric")
  {
    if(nq==0)
    {
      for(i in 1:n)
      {
        AA[i,i:n] = t(B[1:(n-i+1)])
        if(root)
          AA[i,n] = R[n+1-i]/(2*pi)
      }
      AA[1,1] = AA[n,n]
      ##  Use symmetry to construct bottom 'half' of AA
      AAu = AA
      AAu[!upper.tri(AAu)] <- 0
      AA = AA + flipud(fliplr(AAu))
    }
    else
    {
      ## CONSTRUCT THE A MATRIX size n x n
      A = Abuild(n,nq,g,root)
      Ainv = solve(A)
      ## CONSTRUCT THE d MATRIX size n x 1
      for(np in 0:ceiling(n/2-1))
      {
        d = matrix(0,n,1)
        ii = 0
        for(jj in (np-root):(np+1+root-n))
        {
          ii = ii+1
          d[ii] = Bge(jj,nq,B,cc)
        }
        if (root == 1)
          d[n-1] = R[n-np]
        ##  COMPUTE Bhat = inv(A)*d
        Bhat = Ainv%*%d
        AA[np+1,] = t(Bhat)
      }
      ##  Use symmetry to construct bottom 'half' of AA
      AA[(ceiling(n/2)+1):n,] = flipud(fliplr(AA[1:floor(n/2),]))
    }
  }


  ###  symmetric filter

  if (type == "symmetric")
  {
    if(nq==0)
    {
      for(i in 2:ceiling(n/2))
      {
        np = i-1
        AA[i,i:(i+np)] = t(B[1:(1+np)])
        if(root)
          AA[i,i+np] = R[np+1]/(2*pi);
        AA[i,(i-1):(i-np)] = AA[i,(i+1):(i+np)];
      }
      ##  Use symmetry to construct bottom 'half' of AA
      AA[(ceiling(n/2)+1):n,] = flipud(fliplr(AA[1:floor(n/2),]))
    }
    else
    {
      for(np in nq:ceiling(n/2-1))
      {
        nf = np
        nn = 2*np+1
        ## CONSTRUCT THE A MATRIX size nn x nn
        A = Abuild(nn,nq,g,root)
        Ainv = solve(A)
        ## CONSTRUCT THE d MATRIX size nn x 1
        d = matrix(0,nn,1)
        ii = 0
        for(jj in (np-root):(-nf+root))
        {
          ii = ii+1
          d[ii] = Bge(jj,nq,B,cc)
        }
        if(root)
          d[nn-1] = R[nf+1]
        ##  COMPUTE Bhat = inv(A)*d
        Bhat = Ainv%*%d
        AA[np+1,1:(2*np+1)] = t(Bhat)
      }
      ##  Use symmetry to construct bottom 'half' of AA
      AA[(ceiling(n/2)+1):n,] = flipud(fliplr(AA[1:floor(n/2),]))
    }
  }


  ###  fixed length symmetric filter

  if (type == "fixed")
  {
    if(nq==0)
    {
      bb = matrix(0,2*nfix+1,1)
      bb[(nfix+1):(2*nfix+1)] = B[1:(nfix+1)]
      bb[nfix:1] = B[2:(nfix+1)]
      if(root)
      {
        bb[2*nfix+1] = R[nfix+1]/(2*pi)
        bb[1] = R[nfix+1]/(2*pi)
      }
      for(i in (nfix+1):(n-nfix))
        AA[i,(i-nfix):(i+nfix)] = t(bb)
    }
    else
    {
      nn = 2*nfix+1
      ## CONSTRUCT THE A MATRIX size nn x nn
      A = Abuild(nn,nq,g,root)
      Ainv = solve(A)
      ## CONSTRUCT THE d MATRIX size nn x 1
      d = matrix(0,nn,1)
      ii = 0
      for(jj in (nfix-root):(-nfix+root))
      {
        ii = ii+1
        d[ii] = Bge(jj,nq,B,cc)
      }
      if(root)
        d[nn-1] = R[nn-nfix]
      ##  COMPUTE Bhat = inv(A)*d
      Bhat = Ainv%*%d
      for(ii in (nfix+1):(n-nfix))
        AA[ii,(ii-nfix):(ii+nfix)] = t(Bhat)
    }
  }

  ###  Baxter-King filter

  if (type == "baxter-king")
  {
    bb = matrix(0,2*nfix+1,1)
    bb[(nfix+1):(2*nfix+1)] = B[1:(nfix+1)]
    bb[nfix:1] = B[2:(nfix+1)]
    bb = bb - sum(bb)/(2*nfix+1)
    for(i in (nfix+1):(n-nfix))
      AA[i,(i-nfix):(i+nfix)] = t(bb)
  }

  ###  Trigonometric Regression filter

  if(type == "trigonometric")
  {
    jj = 1:(n/2)
    ## find frequencies in desired band omitting n/2;
    jj = jj[((n/pu)<=jj & jj<=(n/pl) & jj<(n/2))]
    if(!any(jj))
      stop("frequency band is empty in trigonometric regression")
    om = 2*pi*jj/n
    if(pl > 2)
    {
      for(t in 1:n)
      {
        for(k in n:1)
        {
          l = t-k
          tmp = sum(cos(om*l))
          AA[t,k] = tmp
        }
      }
    }
    else
    {
      for(t in 1:n)
      {
        for(k in n:1)
        {
          l = t-k
          tmp = sum(cos(om*l))
          tmp2 = (cos(pi*(t-l))*cos(pi*t))/2
          AA[t,k] = tmp + tmp2
        }
      }
    }
    AA = AA*2/n
  }


  ###  check that sum of all filters equal 0 if assuming unit root

  if(root)
  {
    tst = max(abs(c(apply(AA,1,sum))))
    if((tst > 1.e-09) && root)
    {
      warning("Bhat does not sum to 0 ")
      cat("test =",tst,"\n")
    }
  }

  ###	compute filtered time series using selected filter matrix AA

  if(drift)
    x = undrift(x)

  x.cycle = AA%*%as.matrix(x)

  if(type=="fixed" || type=="symmetric" || type=="baxter-king")
  {
    if(nfix>0)
      x.cycle[c(1:nfix,(n-nfix+1):n)] = NA
  }
  x.trend = x-x.cycle
  if(is.ts(xo))
  {
    tsp.x = tsp(xo)
    x.cycle=ts(x.cycle,start=tsp.x[1],frequency=tsp.x[3])
    x.trend=ts(x.trend,start=tsp.x[1],frequency=tsp.x[3])
    x=ts(x,start=tsp.x[1],frequency=tsp.x[3])
  }

  if(type=="asymmetric")
    title = "Chiristiano-Fitzgerald Asymmetric Filter"
  if(type=="symmetric")
    title = "Chiristiano-Fitzgerald Symmetric Filter"
  if(type=="fixed")
    title = "Chiristiano-Fitzgerald Fixed Length Filter"
  if(type=="baxter-king")
    title = "Baxter-King Fixed Length Filter"
  if(type=="trigonometric")
    title = "Trigonometric Regression Filter"

  res <- list(cycle=x.cycle,trend=x.trend,fmatrix=AA,title=title,
              xname=xname,call=as.call(match.call()),
              type=type,pl=pl,pu=pu,nfix=nfix,root=root,drift=drift,
              theta=theta,method="cffilter",x=x)

  return(structure(res,class="mFilter"))
}


###======================================================================
###				Functions
###======================================================================

Bge <- function(jj,nq,B,cc)
{
  ###
  ###  closed form solution for integral of B(z)g(z)(1/z)^j  (eqn 16)
  ###     nq > 0, jj >= 0
  ###     if nq = 0, y = 2*pi*B(absj+1)*cc(1);
  ###
  absj =abs(jj)
  if(absj >= nq)
  {
    dj = B[absj+1]*cc[1] + t(B[(absj+2):(absj+nq+1)])%*%cc[2:(nq+1)]
    dj = dj + t(flipud(B[(absj-nq+1):absj]))%*%cc[2:(nq+1)]
  }
  else if(absj >= 1)
  {
    dj = B[absj+1]*cc[1] + t(B[(absj+2):(absj+nq+1)])%*%cc[2:(nq+1)]
    dj = dj + t(flipud(B[1:absj]))%*%cc[2:(absj+1)]
    dj = dj + t(B[2:(nq-absj+1)])%*%cc[(absj+2):(nq+1)]
  }
  else
    dj = B[absj+1]*cc[1] + 2*t(B[2:(nq+1)])%*%cc[2:(nq+1)]

  y = 2*pi*dj
  return(y)
}

###
###-----------------------------------------------------------------------
###

Abuild <- function(nn,nq,g,root)
{
  ###
  ###  builds the nn x nn A matrix in A.12
  ###   if root == 1 (unit root)
  ###      Abig is used to construct all but the last 2 rows of the A matrix
  ###   elseif root == 0 (no unit root)
  ###      Abig is used to construct the entire A matrix
  ###
  if(root)
  {
    Abig=matrix(0,nn,nn+2*(nq-1))
    for(j in 1:(nn-2))
      Abig[j,j:(j+2*nq)] = t(g)
    A = Abig[,nq:(nn+nq-1)]
    ##   construct A(-f)
    Q = -matrix(1,nn-1,nn)
    ##Q = tril(Q);
    Q[upper.tri(Q)] <- 0
    Fm = matrix(0,1,nn-1)
    Fm[(nn-1-nq):(nn-1)] = g[1:(nq+1)]
    A[(nn-1),] = Fm%*%Q
    ##    construct last row of A
    A[nn,] = matrix(1,1,nn)
  }
  else
  {
    Abig=matrix(0,nn,nn+2*(nq-0))
    for(j in 1:nn)
    {
      Abig[j,j:(j+2*nq)] = c(g)
    }
    A = Abig[,(nq+1):(nn+nq)]
  }
  ##    multiply A by 2*pi
  A = 2*pi*A
  return(A)
}
###
###-----------------------------------------------------------------------
###

undrift <- function(x)
{
  ###
  ###  This function removes the drift or a linear time trend from a time series using the formula
  ###			drift = (x(n) - x(1)) / (n-1).
  ###
  ###  Input:  x - data matrix x where columns represent different variables, x is (n x # variables).
  ###  Output: xun - data matrix same size as x with a different drift/trend removed from each variable.
  ###
  x = as.matrix(x)
  nv = dim(x)
  n = nv[1]
  nvars = nv[2]
  xun = matrix(0,n,nvars)
  dd = as.matrix(0:(n-1))
  for(ivar in 1:nvars)
  {
    drift = (x[n,ivar]-x[1,ivar]) / (n-1)
    xun[,ivar] = x[,ivar] - dd*drift
  }
  if(nvars==1) xun = c(xun)

  return(xun)
}

###
###-----------------------------------------------------------------------
###

###function that reverses the columns of a matrix (matlab equivalent)
flipud <- function(x) {apply(as.matrix(x),2,rev)}

###function that reverses the rows of a matrix (matlab equivalent)
fliplr <- function(x) {t(apply(as.matrix(x),1,rev))}



# A flowFrame object is splitted in bins that are indicate the flow
# rate over time.
# @param second_fraction # change the fraction of seconds according
# to your desire
# @return The returned value is a list with the following components.
# @return \item{anoms}{Data frame containing index, values, and
# expected values.}
#
flow_rate_bin <- function(x, second_fraction = 0.1, timeCh = timeCh,
                          timestep = timestep){

  xx <- exprs(x)[, timeCh]
  idx <- c(1:nrow(x))

  endsec <- ceiling(timestep * max(xx))  # total seconds of the experiment

  lenx <- length(xx)  # num of time ticks

  tbins <- seq(0, endsec/timestep, by = as.numeric(second_fraction)/timestep)  # time bins
  secbin <- seq(0, endsec, by = as.numeric(second_fraction))  # bin expressed in seconds
  minbin <- round(secbin/60, 3)  # bin expressed in minutes
  nrBins <- length(tbins) - 1
  tbCounts <- c(0, hist(xx, tbins, plot = FALSE)$counts)  # number of events per time bin
  expEv <- lenx/(nrBins)  ##median(tbCounts) # expected number of events per bin
  binID <- do.call(c, mapply(rep, x = 1:length(tbCounts), times = tbCounts,
                             SIMPLIFY = FALSE))

  if (length(idx) != length(binID))
    stop("length of cell ID not equal length of bin ID")

  timeFlowData <- list(frequencies = cbind(tbins, minbin, secbin, tbCounts),
                       cellBinID = data.frame(cellID = idx, binID = binID),
                       info = data.frame(second_fraction = second_fraction,
                                         expFrequency = expEv, bins = nrBins))
  return(timeFlowData)
}


# Detection of anomalies in the flow rate using the algorithm
# implemented in the package AnomalyDetection.
flow_rate_check_a <- function(x, FlowRateData, alpha = alpha, use_decomp = use_decomp) {

  fr_frequences <- FlowRateData$frequencies
  fr_cellBinID <- FlowRateData$cellBinID
  second_fraction <- FlowRateData$info["second_fraction"]


  if (length(unique(fr_frequences[, 2])) == 1) {
    fr_autoqc <- NULL
  } else {
    fr_autoqc <- anomaly_detection(fr_frequences[, "tbCounts"], alpha = alpha, use_decomp = use_decomp)
  }

  if (is.null(fr_autoqc) || is.null(fr_autoqc$anoms)) {
    badPerc <- 0
    newx <- x
    goodCellIDs <- fr_cellBinID$cellID
    badCellIDs <- NULL
  } else {
    goodCellIDs <- fr_cellBinID$cellID[!(fr_cellBinID$binID %in% fr_autoqc$anoms$index)]
    badCellIDs <- setdiff(fr_cellBinID$cellID, goodCellIDs)
    badPerc <- round(1 - (length(goodCellIDs)/nrow(fr_cellBinID)), 4)
    params <- parameters(x)
    keyval <- keyword(x)
    sub_exprs <- exprs(x)
    sub_exprs <- sub_exprs[goodCellIDs, ]
    newx <- flowFrame(exprs = sub_exprs, parameters = params, description = keyval)
  }
  cat(paste0(100 * badPerc, "% of anomalous cells detected in the flow rate check. \n"))
  return(list(anoms = fr_autoqc$anoms, frequencies = fr_frequences,
              FRnewFCS = newx,
              goodCellIDs = goodCellIDs, badCellIDs = badCellIDs,
              res_fr_QC = data.frame(second_fraction = second_fraction,
                                     num_obs = fr_autoqc$num_obs, badPerc = badPerc)))
}


# Plot frequency values for a list y, containing the outputs from
# the function flow_rate_check
flow_rate_plot <- function(FlowRateQC) {

  second_fraction <- FlowRateQC$res_fr_QC$second_fraction
  num_obs = FlowRateQC$res_fr_QC$num_obs
  frequencies = as.data.frame(FlowRateQC$frequencies)
  anoms = as.data.frame(FlowRateQC$anoms)
  anoms_points = as.data.frame(cbind(sec_anom = frequencies$secbin[anoms$index], count_anom = anoms$anoms))


  xgraph <- ggplot(frequencies, aes_string(x="secbin", y="tbCounts")) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       text=element_text(size = 14)) + geom_line(colour = "red" )
  xgraph <- xgraph + labs(x= "Seconds", y= paste0("Number of events per 1/",
                                                  1  /second_fraction, " of a second"), title= "Flow Rate")

  # Add anoms to the plot as circles.
  if(!is.null(anoms_points)){
    xgraph <- xgraph + geom_point(data=anoms_points, aes_string(x= "sec_anom", y= "count_anom"), color = "green4", size = 3, shape = 1)
  }

  return(xgraph)
}

# Retrieves the number of events for each FCS file and creates
# a bar plot.
#
flow_set_qc <- function(set){
  if (!is(set, "flowSet"))
    stop("'set' needs to be of class 'flowSet'")

  cells <- as.numeric(fsApply(set, nrow))
  cells[cells == 1] <- NA
  samples <- factor(flowCore::sampleNames(set), levels = flowCore::sampleNames(set))
  cnframe <- data.frame(sampleName = samples , cellNumber = cells)
  return(cnframe)
}

# produce a bar plot
flow_set_plot <- function(N_cell_set, area){

  ggplot(N_cell_set, aes_string(x="sampleName", y = "cellNumber")) +
    geom_bar(stat = "identity",fill= area) + theme_classic() +
    theme( legend.position="none", axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# The events on the upper margin and the outlier in the negative
# range of values are detected and removed.
#
flow_margin_check_a <- function(x,  margin_channels = NULL,
                              side = "both", neg_values = 1) {

  if (is.null(margin_channels)) {
    teCh <- grep("Time|time|TIME|Event|event|EVENT", colnames(x), value = TRUE)
    parms <- setdiff(colnames(x), teCh)
  } else {
    if (!all(margin_channels %in% colnames(x)))
      stop("Invalid channel(s)")
    parms <- margin_channels
  }
  scatter_parms <- grep("FSC|SSC", parms, value = TRUE)

  xx <- c(1:nrow(x))
  yy <- x@exprs[, parms]
  range <- range(x)
  lenx <- length(xx)

  ## lower check
  if ((side == "lower" || side == "both") && neg_values == 1) {
    out_neg_range <- apply(yy, 2, function(x) {
      neg <- which(x < 0)
      # Zscores <- (0.6745*(x[neg] + median(x[neg])))/mad(x[neg]) ## it
      # calculates the Zscore outneg <- neg[which(Zscores < -3.5)]
      min_value <- -3.5 * mad(x[neg]) / 0.6745 + median(x[neg])  # -3.5 is the default threshold
      if (is.na(min_value)) {
        min(x) - 1
      } else {
        min_value
      }
    })
  }

  # n. bad cells for each channel
  if ((side == "lower" || side == "both") && neg_values == 1) {
    neg_bad_len <- sapply(parms, function(x) length(xx[yy[, x] <= out_neg_range[x]]))
  }
  if ((side == "lower" || side == "both") && neg_values == 2) {
    neg_bad_len <- sapply(parms, function(x) length(xx[yy[, x] <= range[1, x]]))
  }
  if (side == "upper" || side == "both") {
    pos_bad_len <- sapply(parms, function(x) length(xx[yy[, x] >= range[2, x]]))
  }

  # badcellIDs
  if ((side == "lower" || side == "both") && neg_values == 1) {
    lowID <- do.call(c, lapply(parms, function(ch) {
      xx[yy[, ch] <= out_neg_range[ch]]
    }))
    if(length(scatter_parms) != 0){   ### check for values less than 0 in scatter parameters
      minSc <- apply(yy[,scatter_parms], 1, function(x){
        min(x)
      })
      low_scatter_ID <- which(minSc < 0)
      lowID <- c(lowID, low_scatter_ID)
    }
  }
  if ((side == "lower" || side == "both") && neg_values == 2) {
    lowID <- do.call(c, lapply(parms, function(ch) {
      xx[yy[, ch] <= range[1, ch]]
    }))
  }
  if (side == "upper" || side == "both") {
    upID <- do.call(c, lapply(parms, function(ch) {
      xx[yy[, ch] >= range[2, ch]]
    }))
  }

  if (side == "lower") {
    summary_bad_cells <- data.frame(lower_range = c(neg_bad_len,
                                                    total_SUM = length(lowID), total_UNIQUE = length(unique(lowID))))
    bad_lowerIDs <- unique(lowID)
    bad_upperIDs <- NULL
    badCellIDs <- unique(lowID)
  } else if (side == "upper") {
    summary_bad_cells <- data.frame(upper_range = c(pos_bad_len,
                                                    total_SUM = length(upID), total_UNIQUE = length(unique(upID))))
    bad_lowerIDs <- NULL
    bad_upperIDs <- unique(upID)
    badCellIDs <- unique(upID)
  } else {
    summary_bad_cells <- data.frame(lower_range = c(neg_bad_len,
                                                    total_SUM = length(lowID), total_UNIQUE = length(unique(lowID))),
                                    upper_range = c(pos_bad_len,
                                                    total_SUM = length(upID), total_UNIQUE = length(unique(upID))))
    bad_lowerIDs <- unique(lowID)
    bad_upperIDs <- unique(upID)
    badCellIDs <- unique(c(lowID,upID))
  }

  goodCellIDs <- setdiff(xx, badCellIDs)
  badPerc <- round(length(badCellIDs)/lenx, 4)

  cat(paste0(100 * badPerc, "% of anomalous cells detected in the dynamic range check. \n"))

  params <- parameters(x)
  keyval <- keyword(x)
  sub_exprs <- exprs(x)
  sub_exprs <- sub_exprs[goodCellIDs, ]
  newx <- flowFrame(exprs = sub_exprs, parameters = params,
                    description = keyval)

  return(list(FMnewFCS = newx, goodCellIDs = goodCellIDs,
              bad_lowerIDs = bad_lowerIDs, bad_upperIDs = bad_upperIDs,
              margin_events = summary_bad_cells, badPerc = badPerc,
              events = lenx))
}


###  graph showing where the anomalies mostly happened
flow_margin_plot <- function(FlowMarginData, binSize = 500) {

  tot_events <- FlowMarginData$events
  bad_lowerIDs <- FlowMarginData$bad_lowerIDs
  bad_upperIDs <- FlowMarginData$bad_upperIDs

  if (missing(binSize) || is.null(binSize) || is.na(binSize))
    binSize <- 500
  nrBins <- floor(tot_events/binSize)

  cf <- c(rep(1:nrBins, each = binSize), rep(nrBins + 1, tot_events - nrBins * binSize))
  tmpx <- split(1:tot_events, cf)

  if(length(bad_lowerIDs) != 0 && length(bad_upperIDs) != 0){
    lowline <- sapply(tmpx, function(x){
      length(which(bad_lowerIDs %in% x))
    })
    upline <- sapply(tmpx, function(x){
      length(which(bad_upperIDs %in% x))
    })
    ymax <- max(lowline, upline)
    plot(lowline, type ="l", col = "blue", bty ="n",
         ylim = c(0, ymax), xlab = "Bin ID",
         ylab = "Number of events removed", cex.lab=1 )
    lines(upline, col = "red")
    legend("top", c("Negative Outliers", "Upper Margin Events"), lty = 1,bty = "n", cex = 1,
           col = c("blue", "red"))
  }else if( length(bad_lowerIDs) != 0 && length(bad_upperIDs) == 0){
    lowline <- sapply(tmpx, function(x){
      length(which(bad_lowerIDs %in% x))
    })
    plot(lowline, type ="l", col = "blue", bty ="n", xlab = "Bin ID",
         ylab = "Number of events removed", cex.lab=1 )
    legend("top", c("Negative Outliers"), lty = 1,bty = "n", cex = 1,
           col = "blue")
  }else if( length(bad_lowerIDs) == 0 && length(bad_upperIDs) != 0){
    upline <- sapply(tmpx, function(x){
      length(which(bad_upperIDs %in% x))
    })
    plot(upline, type ="l", col = "red", bty ="n", xlab = "Bin ID",
         ylab = "Number of events removed", cex.lab=1 )
    legend("top", c("Upper Margin Events"), lty = 1,bty = "n", cex = 1,
           col = "red")
  }
}

# A flowFrame object is splitted in bins with equal number of events
# and for each bin the median is calculated.
#
# @param x: a flowFrame object
# @param channels: channel names for the selected markers
# @param binSize: the size of bin
flow_signal_bin_a <- function(x, channels = NULL, binSize = 500,
                            timeCh = timeCh, timestep = timestep, TimeChCheck = TimeChCheck) {

  ## some sanity checking
  if (!is(x, "flowFrame"))
    stop("'x' needs to be of class 'flowFrame'")

  if (is.null(channels) || missing(channels) || is.na(channels)) {
    parms <- setdiff(colnames(x), timeCh)
  } else {
    if (!all(channels %in% colnames(x)))
      stop("Invalid channel(s)")
    parms <- channels
  }

  ### Retriving time and expression info
  exp <- exprs(x)
  if (!is.null(TimeChCheck)) {
    timex <- seq(from = 0, length.out = nrow(x), by = 0.1)
  }else{
    timex <- exp[, timeCh]
  }
  yy <- exp[, parms]  # channels data
  idx <- c(1:nrow(x))
  seconds <- timex * timestep
  lenSec <- length(seconds)  # num of time ticks
  uniSeconds <- unique(seconds)  # num of unique time tick
  nrBins <- floor(lenSec/binSize)  # num of bins

  cf <- c(rep(1:nrBins, each = binSize), rep(nrBins + 1, lenSec - nrBins * binSize))  # id bins
  stopifnot(length(cf) == lenSec)
  tmpx <- split(seconds, cf)
  xx2 <- sapply(tmpx, mean)  # mean of each time bin  (x axis)
  yy2 <- as.matrix(ddply(as.data.frame(yy), .(cf), colwise(median)))[, -1]

  return(list(exprsBin = cbind(timeSec = xx2, yy2), cellBinID = data.frame(cellID = idx, binID = cf),
              bins = length(unique(cf)), binSize = binSize))
}

#########################################################################
# Detection of shifts in the median intensity signal detected
# by the laser of the flow cytometry over time
flow_signal_check_a <- function(x, FlowSignalData, ChannelRemove = NULL,
                              pen_valueFS = pen_valueFS, maxSegmentFS = maxSegmentFS, outlier_remove = FALSE) {

  fs_cellBinID <- FlowSignalData$cellBinID
  fs_res <- FlowSignalData$exprsBin
  ### log transformation.
  # fs_res[which(fs_res <= 1 & fs_res >= -1)] <- 0
  # fs_res[which(fs_res > 1)] <- log(fs_res[which(fs_res > 1)])
  # fs_res[which(fs_res < -1)] <- -log(abs(fs_res[which(fs_res < -1)]))

  teCh <- grep("Time|time|TIME|Event|event|EVENT", colnames(fs_res), value = TRUE)
  parms <- setdiff(colnames(fs_res), teCh)

  ##### scale and sum the value of each channel and find the outliers
  scale_sign <- apply(fs_res[, parms],2, scale)
  sum_sign <- apply(abs(scale_sign),1, sum)
  outup_tr <- (+3.5 * mad(sum_sign) + (0.6745 * median(sum_sign)))/0.6745
  FS_out <- which(sum_sign >outup_tr)

  #### Remove channel from the changepoint analysis
  if (!is.null(ChannelRemove)) {
    ChannelRemove_COMP <- grep(paste(ChannelRemove, collapse="|"),
                               colnames(fs_res), value = TRUE)
    #  cat(paste0("The channels whose signal acquisition will not be checked are: ",
    #   paste(ChannelRemove_COMP, collapse = ", "), ". \n"))
    parms <- setdiff(parms, ChannelRemove_COMP)
  }

  if(outlier_remove){
    ##transform outliers to median values
    if( length(FS_out) != 0 ){
      fs_res_adj <- apply(fs_res[, parms], 2, function(x) {
        med <- median(x)
        x[FS_out] <- med
        return(x)
      })
      badPerc_out <- round((length(FS_out)/nrow(fs_res)),4)
      # cat(paste0(badPerc_out * 100, "% of outliers found in channels' signal. \n"))
    }else{
      fs_res_adj <- fs_res[,parms]
      #  cat("0% of outliers found in channels' signal. \n")
      badPerc_out <- 0
    }

    cpt_res <- suppressWarnings(cpt.meanvar(t(fs_res_adj),
                                            pen.value = pen_valueFS, Q = maxSegmentFS,
                                            penalty =  "Manual" , test.stat = "Normal",
                                            method = "BinSeg", param.estimates = FALSE))
  }else{

    cpt_res <- suppressWarnings(cpt.meanvar(t(fs_res[, parms]),
                                            pen.value = pen_valueFS, Q = maxSegmentFS,
                                            penalty =  "Manual" , test.stat = "Normal",
                                            method = "BinSeg", param.estimates = FALSE))
    badPerc_out <- 0
  }

  list_seg <- lapply(1:length(cpt_res), function(x) {
    # it retrieves the changepoints that has been detected for each channel
    cpt_res[[x]]@cpts
  })
  names(list_seg) <- parms
  list_seg <- as.list(list_seg)

  list_cpt <- union(1, sort(unique(unlist(list_seg[parms]))))  # retrieve and order all the changepoints
  diff_cpt <- sapply(2:length(list_cpt), function(n) {
    # calculate the difference between each segment
    list_cpt[n] - list_cpt[n - 1]
  })
  max_diff <- which(diff_cpt == max(diff_cpt))
  max_seg <- c(list_cpt[max_diff], list_cpt[max_diff + 1] )  # selecting the biggest segment

  list_seg <- lapply(list_seg, function(x) setdiff(x, nrow(fs_res)))

  len_cpt <- sapply(list_seg, length)
  nam_cpt <- gsub("<|>","",names(len_cpt))
  nozero_cpt <- as.numeric(which(len_cpt != 0))
  zero_cpt <- as.numeric(which(len_cpt == 0))
  if(length(nozero_cpt) == 0){
    ch_no_cpt <- nam_cpt[zero_cpt]
    tab_cpt <- NULL
  }else{
    # cat(paste("Changepoint(s) detected in the channels: ",
    # paste(names(len_cpt[nozero_cpt]), collapse = ", "), sep = ""), fill = TRUE)
    ch_cpt <- list_seg[nozero_cpt]
    ch_no_cpt <- nam_cpt[zero_cpt]

    max_n_cpt <- max(sapply(ch_cpt, length))
    tab_cpt <- ldply(ch_cpt, function(x) c(x, rep(NA, max_n_cpt - length(x))),
                     .id = NULL)
    rownames(tab_cpt) <- nam_cpt[nozero_cpt]
    tab_cpt <- as.matrix(tab_cpt)
    tab_cpt[which(is.na(tab_cpt))] <- ""
    colnames(tab_cpt) <- 1:length(tab_cpt[1, ])
  }
  # percentage bad cell detected with the changepoint method
  badPerc_cp <- round(1 - ((max_seg[2] - max_seg[1])/(length(fs_res[, 1]) - 1)),4)

  cat(paste0(100 * badPerc_cp, "% of anomalous cells detected in signal acquisition check. \n"))

  # retrieve ID of good cells
  if(outlier_remove){
    fs_cellBinID <- fs_cellBinID[which(!fs_cellBinID[, 2] %in% FS_out),]
  }
  goodCellIDs <- fs_cellBinID[which(fs_cellBinID[, 2] >= max_seg[1] &
                                      fs_cellBinID[, 2] <= max_seg[2]), 1]

  badPerc_tot <- round(1 - length(goodCellIDs)/nrow(fs_cellBinID),4)

  params <- parameters(x)
  keyval <- keyword(x)
  sub_exprs <- exprs(x)
  sub_exprs <- sub_exprs[goodCellIDs, ]  ## check if the Id Correspond!
  newx <- flowFrame(exprs = sub_exprs, parameters = params, description = keyval)

  return(list(FSnewFCS = newx, exprsBin = FlowSignalData$exprsBin, Perc_bad_cells = data.frame(badPerc_tot,badPerc_cp, badPerc_out),
              goodCellIDs = goodCellIDs, tab_cpt = tab_cpt, ch_no_cpt =ch_no_cpt,
              segm = max_seg, FS_out = FS_out, outlier_remove = outlier_remove))
}


# Plot the flourescence intensity for each channel of a flowFrame
# over time, highlighting the wider segment that do not show shifts
# of the median intensity
# @param exprsBin give the exprsBin object from FlowSignalData
# @param segm give the segm object from FlowSignalQC
# @param FS_out give the FS_out object from FlowSignalQC
flow_signal_plot <- function(FlowSignalQC) {

  exprsBin <- FlowSignalQC$exprsBin
  segm <- FlowSignalQC$segm
  FS_out <- FlowSignalQC$FS_out
  outlier_remove <- FlowSignalQC$outlier_remove

  binID <- 1:nrow(exprsBin)
  teCh <- grep("Time|time|Event|event", colnames(exprsBin), value = TRUE)
  parms <- setdiff(colnames(exprsBin), teCh)
  dataORIG <- exprsBin[, parms]  # first channel is time

  if(length(FS_out) != 0){
    data <- apply(dataORIG, 2, function(x){
      overMAX <- FS_out[x[FS_out] > max(x[-FS_out])]
      overMIN <- FS_out[x[FS_out] < min(x[-FS_out])]
      x[overMAX] <- max(x[-FS_out])
      x[overMIN] <- min(x[-FS_out])
      return(x)
    })
    data <- as.data.frame(data)
  }else{
    data <- as.data.frame(dataORIG)
  }
  data$binID <- binID
  longdata <- melt(data, id.vars = "binID", variable.name = "marker",
                   value.name = "value")
  FS_graph <- ggplot(longdata, aes_string(x = "binID", y = "value", col = "marker"),
                     environment = environment()) + labs(x = "Bin ID",
                                                         y = "Median Intensity value") +
    facet_grid(marker ~ ., scales = "free") + theme_bw() +
    theme(strip.text.y = element_text(angle = 0,
                                      hjust = 1), axis.text.y = element_text(size = 6),
          legend.position = "none") +
    scale_x_continuous(breaks= pretty_breaks(n =10)) +
    scale_y_continuous(breaks= pretty_breaks(n =3)) +
    geom_rect(aes(xmin = segm[1], xmax = segm[2], ymin = -Inf,
                  ymax = Inf), fill = "khaki1", linetype = 0) + geom_line()
  # Add anoms to the plot as circles.
  if(outlier_remove){
    longdata_out <- melt(data[FS_out,], id.vars = "binID",
                         variable.name = "marker",value.name = "value")
    FS_graph <- FS_graph + geom_point(data=longdata_out, aes_string(x= "binID",
                                                                    y= "value"), color = "green4", size = 2, shape = 1)
  }
  return(FS_graph)
}

## remove last slash if present
strip.sep <- function(name) {
  ifelse(substr(name,nchar(name),nchar(name))==.Platform$file,
         substr(name,1,nchar(name)-1),name)
}

# Guess which channel captures time in a exprs, flowFrame or flowset
findTimeChannel <- function(xx) {
  time <- grep("^Time$", colnames(xx), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(time)) {
    if (is(xx, "flowSet") || is(xx, "ncdfFlowList"))
      xx <- exprs(xx[[1]]) else if (is(xx, "flowFrame"))
        xx <- exprs(xx)
      cont <- apply(xx, 2, function(y) all(sign(diff(y)) >= 0))
      time <- names(which(cont))
  }
  if (!length(time) || length(time) > 1)
    time <- NULL
  return(time)
}

# Check if the Fcs file is ordered according to time otherwise it order it.
ord_fcs_time <- function(x, timeCh){

  xord <- order(exprs(x)[, timeCh])

  if( !identical(xord, 1:nrow(x)) ){
    warning(paste0("Expression data in the file ", basename(keyword(x)$FILENAME),
                   " were not originally ordered by time."))
    params <- parameters(x)
    keyval <- keyword(x)
    sub_exprs <- exprs(x)[xord, ]
    newx <- flowFrame(exprs = sub_exprs, parameters = params,
                      description = keyval)
    return(newx)
  }else{
    return(x)
  }
}

## create new flowFrame with the parameter indicating good and bad cells
addQC_a <- function(QCvector, remove_from, sub_exprs, params, keyval){

  rs <- attr(sub_exprs, "ranges")
  rs <- c(rs, rs[1])
  sub_exprs <- cbind(sub_exprs, QCvector)
  attr(sub_exprs, "ranges") <- rs
  NN <- as.numeric(keyval["$PAR"]) + 1
  names(dimnames(sub_exprs)[[2]]) <- sprintf("$P%sN", 1:NN)
  pnr <- paste0("$P", NN, "R")
  pnb <- paste0("$P", NN, "B")
  pne <- paste0("$P", NN, "E")
  pnn <- paste0("$P", NN, "N")
  pns <- paste0("$P", NN, "S")
  flowCorePnRmax <- paste0("flowCore_$P", NN, "Rmax")
  flowCorePnRmin <- paste0("flowCore_$P", NN, "Rmin")
  o <- params@data
  o[length(o[,1]) + 1,] <- c(paste0("remove_from_", remove_from), "QC", as.numeric(keyval$`$P1R`), 0, 20000)
  rownames(o)[length(o[,1])] <- paste("$P", NN, sep = "")

  outFCS <- new("flowFrame", exprs=sub_exprs, parameters=new("AnnotatedDataFrame",o), description=keyval)
  description(outFCS)[pnr] <- max(20000, description(outFCS)$`$P1R`)
  description(outFCS)[pnb] <- description(outFCS)$`$P1B`
  description(outFCS)[pne] <- "0,0"
  description(outFCS)[pnn] <- paste0("remove_from_", remove_from)
  description(outFCS)[pns] <- "QC"
  description(outFCS)$`$PAR` <- NN
  description(outFCS)[flowCorePnRmax] <- 20000
  description(outFCS)[flowCorePnRmin] <- 0
  outFCS
}


