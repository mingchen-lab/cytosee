### This part is from clusterCons ### Because this is no longer offered by cran so I determined to integrated this to my project
###############################################################################
# Alogrithm specification file for clusterCons
#
# Author: Dr. T. Ian Simpson
# Affiliation : University of Edinburgh
# E-mail : ian.simpson@ed.ac.uk
# Version : 1.0
###############################################################################
#The following clustering algorithms from cluster are specified here :-
#		agnes  - agglomerative hierarchical clustering
#		pam    - partitioning around medoids, a more robust version of k-means
#		kmeans - k-means clustering
#		hclust - hierarchical cluster analysis on a set of dissimilarities
#		diana  - divisive hierarchical clustering
#		apcluster - affinity propagation clustering
###############################################################################
#load the cluster package
#' @import cluster

#kmeans
kmeans_clmem <- function(x,clnum,params=list()){
  #remove the diss argument as kmeans doesn't take this argument
  params$diss=NULL;
  #set the MacQueen algorithm as default we have seen random failures with Hartigan-Wong
  params$algorithm='Mac';
  #if object is dist then convert back to matrix
  if(class(x)=='dist'){
    x <- as.matrix(x);
  }
  if(length(params)==0){
    params <- list(centers=clnum);
  }
  else{
    #add or override the centers number in existing params to the current clnum (resolves potential conflict)
    params$centers=clnum;
  }
  params$x <- x;
  clmem <- data.frame(do.call(kmeans,params)$cluster);
  names(clmem) <- c('cm');
  return(clmem);
}
#examples
#1. simple
#kmeans_clmem(cluster_example_data2,clnum=4)
#2. with change of kmeans algorithm, note clnum overides centers
#kmeans_clmem(cluster_example_data2,clnum=6,params=list(algorithm="Lloyd",centers=23))

#hclust
hclust_clmem <- function(x,clnum,params=list()){
  #remove any erroneous diss parameters, hclust only takes 'dist' objects
  params$diss=NULL;
  #if the data is not a dissimilarity already make one. The user should not use hclust without providing the
  #distance matrix as it only takes distance matrices.
  if(class(x)!='dist'){
    dm <- dist(x);
    params$d=dm;
  }
  else{
    dm <- x;
    params$d=dm;
  }
  #perform the clustering
  hc <- do.call(hclust,params)
  #grab the membership object
  clmem <- data.frame(cutree(hc,clnum))
  #tidy the column header
  names(clmem) <-c('cm')
  return(clmem);
}

###############################################################################
# clusterCons - Consensus clustering functions for R
#
# Author: Dr. T. Ian Simpson
# Affiliation : University of Edinburgh
# E-mail : ian.simpson@ed.ac.uk
# Version : 1.0
#
###############################################################################

#read in required libraries


#check data integrity
data_check<-function(x){
  #is it a data.frame
  if(is.data.frame(x)!=1){stop('x is not a data frame')}
  #row.names
  if(length(row.names(x))<1){stop('x has no row names defined')}
  #col.names
  if(length(names(x))<=0){stop('x has no column names defined')}
  #values as numeric (as you cannot have mixed class checking 1 column is sufficient)
  if(!is.numeric(x[,1])){stop('x is not numeric')}
  #check for missing values
  if(sum(is.na(x))!=0){stop('x contains missing values')}
  #check that there are at least two columns and two rows
  if(sum(dim(x)>=2)!=2){stop('the data matrix is not at least 2x2, clustering is pointless')}
  return(TRUE);
}

#extract the data.frame from the ExpressionSet so that it can be used by cluscomp
expSetProcess<-function(x){
  #check that x is an ExpressionSet object, just a catchall in case of direct use
  if(class(x)=='ExpressionSet'){
    #extract the expression data from the ExpressionSet object
    x_data <- x@assayData$exprs;
    #check the rows and columns are uniquely identified
    #check that each of the rownames is unique
    if(
      length(unique(rownames(x_data))) == length(rownames(x_data))
      &&
      #check that each of the colnames is unique
      length(unique(colnames(x_data))) == length(colnames(x_data))
    ){
      #now call the data_check
      x_data <- data.frame(x_data);
      #check the names etc in the conversion here
      if(data_check(x_data)==TRUE){
        return(x_data);
      }
      else{
        stop('The expression set has failed the data check, cannot proceed');
      }
    }
    else{
      stop('row and column names are not unique, cannot process');
    }
  }
  else{
    stop('The data object is not an ExpressionSet');
  }
}

#Basic S4 class to carry the consensus matrices and their full clustering result
#where cm = consensus matrix, rm = reference matrix (full clustering result), a = algorithm and k = cluster number
setClass("consmatrix",representation(cm="matrix",rm="data.frame",a="character",k="numeric",call="call"))
#define a validity checking function
validConsMatrixObject <- function(object){
  if(length(object@cm)>0 & length(object@rm)>0 & length(object@a)>0 & length(object@k)>0) TRUE
  else{paste('One of the required variables has zero length')}
}
setValidity("consmatrix",validConsMatrixObject)

setClass("mergematrix",representation(cm="matrix",k="numeric",a="character"))
validMergeMatrixObject <- function(object){
  if(length(object@cm)>0 & length(object@k)>0 & length(object@a)>0 ) TRUE
  else{paste('One of the required variables has zero length')}
}
setValidity("mergematrix",validMergeMatrixObject)

setClass("memroblist",representation(mrl="data.frame"));
validMemRobListObject <- function(object){
  if(length(object@mrl)>0) TRUE
  else{paste('The memebership robustness list is empty')}
}
setValidity("memroblist",validMemRobListObject)

setClass("memrobmatrix",representation(mrm="matrix"))
validMemRobMatrixObject <- function(object){
  if(length(object@mrm)>0) TRUE
  else{paste('The full membership robustness matrix is empty')}
}
setValidity("memrobmatrix",validMemRobMatrixObject)

#New additions class objects for aucs and deltak plot data, simple inheritance objects
#define a class aucs
setClass("auc",contains=c('data.frame'));
#define a validity checking function
validAUCObject <- function(object){
  if(class(object)!='auc'){stop('The object is not a valid AUC object')}
  else{return(TRUE)};
}
setValidity("auc",validAUCObject)

#define a class dk
setClass("dk",contains=c('data.frame'));
#define a validity checking function
validDkObject <- function(object){
  if(class(object)!='dk'){stop('The object is not a valid dk object')}
  else{return(TRUE)};
}
setValidity("dk",validDkObject)

#TODO probably should put class validation methods in here

#' @import lattice

cluscomp<-function(x,diss=FALSE,algorithms=list('kmeans'),alparams=list(),alweights=list(),clmin=2,clmax=10,prop=0.8,reps=50,merge=0){
  #CHECK INPUTS
  #first check if the input data is an ExpressionSet if so extract the data
  if(class(x)=='ExpressionSet'){x <- expSetProcess(x);};
  #check integrity of data (what about NAs)
  if(data_check(x)!=TRUE){stop('the provided data fails the data integrity check')}
  #check that the algorithms are passed as a list
  if(class(algorithms)!='list'){stop('You must pass the algorithm names in a list')};
  #check that all of the algorithms specified are available
  for(i in 1:length(algorithms)){
    if(existsFunction(algorithms[[i]])!=1){stop('One or more of the specified algorithms does not exist')} else {}
  }
  #alweight checking
  #first check for alweights if NULL then create a list of 1's as long as the number of cms
  if(length(alweights)==0){
    alweights = as.list(rep('1',length(algorithms)))
  }
  #otherwise check that number of weighting elements equals the number of algorithms specified
  else{
    if(length(alweights)!= length(algorithms)){stop('If you sepcify algorithm weightings their number must equal the number of algorithms being used')}
  }
  #clmin,clmax,reps integers
  if(clmin < 2){stop('cannot have less than two clusters')}
  #prop between 0 and 1
  if(prop <= 0 | prop >1){stop('resampling proportion must be greater than 0 and less than 1')}
  #check merge
  if(merge != 0 & merge !=1){stop('merge can only be 0 or 1')}
  #normalise algorithm weightings
  total = sum(as.numeric(alweights));
  norm = as.numeric(alweights)/total;
  #if we get past all this now call the re-sampling
  sample_number = dim(x)[1];

  #list to hold all of the consmatrix objects
  cmlist <- list();

  for(clnum in clmin:clmax)
  {
    print(paste('cluster number k=',clnum),quote=FALSE);

    #if merging create a matrix to hold the merge
    if(merge>0){
      mm <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)))
    }
    #do the re-sample

    #for each algorithm
    for(a in 1:length(algorithms)){
      print(paste('running experiment',a,'algorithm:',algorithms[[a]],sep=' '),quote=FALSE);
      #be aware this connectivity matrix is N^2 so don't try to cluster with 10'000s of genes
      final_conn <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)));
      final_id <- final_conn;
      #matrix to hold the clustering results
      #cl_reps = matrix(0,nrow=dim(x)[1],ncol=reps,dimnames=list(row.names(x),1:reps));
      #get the algorithm to call for cm
      algo_clmem <- paste(algorithms[[a]],'_clmem',sep='');
      current_algo <- get(algo_clmem);
      #check if params have been specified
      if(length(alparams)!=0){
        current_params <- alparams[[a]];
        #check if the diss parameter is set and true, in case user forgets diss flag
        if(!is.null(current_params$diss)){
          diss=TRUE;
        }
      }
      else{
        #here someone has specified that x is a distance matrix, but not included it in the params
        if(diss==TRUE){
          current_params = list(diss=TRUE);
        }
        #otherwise x not a distance matrix params stay empty
        else{
          current_params = list();
        }
      }
      for(i in 1:reps){
        #the generic clustering algorithm calling method
        #checking if data is in fact a distance matrix (data.frame)
        #either by flag diss being set or diss being set in parameters
        if (diss==TRUE){
          #we have already put the object through the data_check
          #we need to check that the row and column names are the same in the distance data.frame
          if(unique(names(x)==row.names(x))==FALSE){stop('the row and column names of the distance matrix are not the same or the distance matrix is not square !')}
          #if OK perform the sample
          else{
            samp_row_names <- sample(row.names(x),as.integer(prop*sample_number));
            samp_x <- x[samp_row_names,samp_row_names];
            #convert to object of class 'dissimilarity'
            samp_x <- as.dist(samp_x);
          }
        }
        else{
          samp_x <- x[sample(row.names(x),as.integer(prop*sample_number)),]; #ORIGINAL
        }
        clmem <- current_algo(samp_x,clnum,params=current_params);
        #cl_reps[row.names(clmem),i] <- clmem$cm # to store implictly the cluster membership results
        #now we can grow the connectivity matrix on the fly by iteration
        #create a current instance of the zeroes indexed array (without pseudocount)
        conn <- matrix(0,nrow=dim(x)[1]+1,ncol=dim(x)[1]+1,dimnames=list(c('x',row.names(x)),c('x',row.names(x))))
        #put in the cluster memberships in the x row,col
        conn[row.names(clmem),1] <- as.matrix(clmem)
        conn[1,row.names(clmem)] <- as.matrix(clmem)
        #now preform the population
        #this needs be re-written to include only i>=j i.e. half the compute (then copy over the opposite half)
        for(i in 2:dim(conn)[1]){
          for(j in 2:dim(conn)[1]){
            if((conn[i,1] == conn[1,j]) & conn[i,1] != 0){
              conn[i,j] <- 1
              conn[j,i] <- 1
            }
            else{
              conn[i,j] <-0
              conn[j,i] <-0
            }
          }
        }
        #cut to the connectivity matrix
        final <- conn[2:dim(conn)[1],2:dim(conn)[1]]
        id = final;
        id[] <- 0;
        id[row.names(clmem),row.names(clmem)] <- 1;
        final_conn = final_conn+final;
        final_id = final_id+id;
      }
      #this is the consensus matrix
      consensus = final_conn/final_id;
      #check for NAs (this is safer than pseudo counting, replace NAs with 0 i.e. they are never drawn together and/or never connected)
      ind <- is.na(consensus)
      consensus[ind] <- 0
      #perform the full clustering as reference
      rm <- current_algo(x,clnum,params=current_params)
      #now create the S4 class
      current_consmatrix <- new('consmatrix',cm=consensus,rm=rm,a=algorithms[[a]],k=clnum,call=match.call());
      #add the consmatrix object to the list
      cmlist[paste('e',a,'_',algorithms[[a]],'_k',clnum,sep='')] <- current_consmatrix
      #now add to the running matrix for merge weighting by the correct algorithm weighting if merge !=0
      if(merge!=0){
        weighted <- current_consmatrix@cm*norm[a]
        mm <- mm+weighted
      }
    }
    if(merge!=0){
      mm <- new('mergematrix',cm=mm,k=clnum,a='merge')
      cmlist[paste('merge_k',clnum,sep='')] <- mm
    }
  }
  return(cmlist)
}

#simple example
#cmr <- cluscomp(testdata,clmin=5,clmax=5,prop=0.8,reps=10)
#more complex example
#alp <- list(method='complete')
#cmr <- cluscomp(testdata,algorithms=list('agnes','pam'),alparams=list(alp,list()),clmin=5,clmax=5,prop=0.8,reps=10) # add params to this
#even more complex with full result retention i.e. both consensus and merge matrices
#cmr <- cluscomp(testdata,algorithms=list('kmeans','pam'),alweights=list(1,0),clmin=4,clmax=6,prop=0.8,reps=3,merge=2)

#function to get the cluster robustness from the consensus matrices
clrob <- function(x,rm=data.frame()){
  if(class(x) == 'consmatrix'){
    cmref <- x@rm
  }
  else{
    if(length(rm)==0){stop('You need to specify a reference matrix for a merge consensus matrix')}
    else{
      cmref <- rm
    }
  }

  cmr <- x@cm

  #get the number of clusters from the reference matrix
  clnum <- length(unique(cmref$cm));
  #set up the matrix to hold the robustness values
  cl_rob = data.frame(matrix(0,clnum,1));
  #label up the matrix
  row.names(cl_rob) <- seq(1:clnum);
  names(cl_rob) <- c('rob');


  #populate the robustness matrix
  for(k in 1:clnum){
    #sub matrix by cluster
    #the as.matrix call here forces to matrix in case any "cluster" only has one member
    sm = as.matrix(cmr[row.names(cmref)[cmref$cm==k],row.names(cmref)[cmref$cm==k]])

    cl_sum <- 0;
    cl_size <- dim(sm)[1];

    for(i in 1:cl_size){
      for(j in 1:cl_size){
        if(i<j){
          cl_sum = cl_sum+sm[i,j];
        }
      }
    }
    #denominator
    meas_sum = 1/(cl_size*(cl_size-1)/2);
    #cluster robustness
    curr_cl_rob = cl_sum*meas_sum;
    if(is.na(curr_cl_rob)){
      cl_rob[k,1] = 0
    }
    else{
      cl_rob[k,1] = curr_cl_rob
    }
  }
  return(cl_rob)
}
##example
#cmr <- cluscomp(testdata,algorithms=list('kmeans','pam'),alweights=list(1,0),clmin=4,clmax=6,prop=0.8,reps=3,merge=2)
##look at the contents of cmr
#summary(cmr)
##retreive the cluster robustness for one of the runs in cmr
#cr <- cl_rob(cmr$kmeans_4)

#function to get the member robustness from the consensus matrices
memrob <- function(x,rm=data.frame()){
  if(class(x) == 'consmatrix'){
    cmref <- x@rm
  }
  else{
    if(length(rm)==0){stop('You need to specify a reference matrix for a merge consensus matrix')}
    else{
      cmref <- rm
    }
  }

  consensus <- x@cm

  mem_rob = matrix(0,dim(consensus)[1],length(summary(as.factor(cmref$cm))),dimnames = list(row.names(consensus),1:length(summary(as.factor(cmref$cm)))))

  for(k in 1:length(summary(as.factor(cmref$cm)))){
    for(i in 1:dim(consensus)[1]){
      Ik = row.names(cmref)[cmref$cm==k] #where k is the cluster number

      ind = Ik[Ik != row.names(consensus)[i]] # exclude the index for i = j if it is there

      sigma = apply(as.matrix(consensus[ind,i]),2,sum) #perform the sigma sum on index

      ei = row.names(consensus)[i] # get the current member we are checking

      Nk = summary(as.factor(cmref$cm))[k] # get the current cluster size


      if(sum(ei == Ik) ==1){		#if ei is a member of Ik
        mik = (1/(Nk-1))*sigma
      }
      else{				#if ei is not a member of Ik
        mik = (1/Nk)*sigma
      }
      mem_rob[i,k] = mik
    }
  }
  #what you might want to do here is have the full object in a slot and output the mem_rob for the ref clustering (which is what you actually want)
  mem_rob_list <- list();
  for(i in 1:dim(mem_rob)[2]){
    cl_mem_rob <- (mem_rob[(cmref==i),i])
    current_list <- data.frame(sort(cl_mem_rob,dec=TRUE))
    names(current_list) <- 'mem_rob'
    current_mem_rob_list <- new('memroblist',mrl=current_list);
    mem_rob_list[paste('cluster',i,sep='')] <- current_mem_rob_list
  }
  mem_rob_list['resultmatrix']<- new('memrobmatrix',mrm=mem_rob);
  mem_rob_list['algo']<- x@a;
  mem_rob_list['type']<- class(x);
  return(mem_rob_list)
}
##example
#cmr <- cluscomp(testdata,algorithms=list('kmeans','pam'),alweights=list(1,0),clmin=4,clmax=6,prop=0.8,reps=3,merge=2)
##look at the contents of cmr
#summary(cmr)
##retreive the member robustness for one of the runs in cmr
#mr <- mem_rob(cmr$kmeans_4)
##look at the contents of the membership robustness object
#summary(mr)
##see the membership robustness of cluster1
#mr$cluster1
##extract the slot object (a data.frame)
#cluster1 <- mr$cluster1@mrl

#AUC - calculate the area under the curve
auc <- function(x){
  n = dim(x)[1]

  m = n*(n-1)/2

  xi <- matrix(0,m,1)

  k=1

  for (i in 1:n){
    for(j in 1:n){
      if(i < j){
        xi[k] <- x[i,j]
        k=k+1
      }
    }
  }

  #sort the xi
  xi <- sort(xi,dec=FALSE)


  #now create the CDF(xi)
  cdfxi = matrix(0,m,1);

  for(i in 1:m){
    #not counting not summing (we are summing the T/F i.e. 1/0)
    cdfxi[i] = sum(xi<=xi[i])/m
  }

  #now create the xi-x(i-1) column
  xidiff = matrix(0,(m-1),1)

  for(i in 2:m){
    xidiff[i-1] = xi[i] - xi[i-1]
  }

  #now put together from 2:l
  test <- data.frame(cbind(xi[2:m],cdfxi[2:m],xidiff))
  #
  names(test)<- c('xi','cdfxi','xidiff')
  #
  test$auc <- test$cdfxi*xidiff

  auc <- sum(test$auc)
  return(auc)
}

#calculate the area under the curves from a consensus clustering result object
aucs <- function(x){
  aucs <- NULL
  k <- NULL
  a<- NULL
  for(i in x){
    a <- c(a,toString(i@a))
    k <- c(k,i@k)
    if(class(i)=='consmatrix'){
      print(paste('calculating AUC for consensus matrix',i@a,i@k,sep=' '),q=FALSE)
    }
    else{
      print(paste('calculating AUC for merge matrix',i@k,sep=' '),q=FALSE)
    }
    aucs <- c(aucs,auc(i@cm))
  }
  a <- data.frame(k=as.factor(k),a=as.factor(a),aucs=aucs)
  a <- new('auc',a);
  return(a)
}

#look at the change in the AUC by cluster number
deltak <- function(x){
  #check input is a valid AUC object
  if(validAUCObject(x)==TRUE){
    deltaks <- data.frame()
    for(i in levels(x$a)){
      current_aucs <- x[x$a==i,]
      deltak <- NULL
      k <- as.numeric(levels(current_aucs$k))
      a <- current_aucs$a
      for(j in 1:length(current_aucs$aucs)){
        if(j==1){
          deltak <- c(deltak,current_aucs$aucs[j])
        }
        else{
          deltak <- c(deltak,((current_aucs$aucs[j]-current_aucs$aucs[j-1])/current_aucs$aucs[j-1]))
        }
      }
      current_deltaks <- data.frame(k=as.factor(k),a=as.factor(a),deltak=deltak)
      deltaks <- rbind(deltaks,current_deltaks);
    }
    deltaks <- new('dk',deltaks);
    return(deltaks)
  }
}




############ this part is for parallel computing for consensus ###########################

#' @import foreach
#' @import doParallel
#' @import iterators
#' @export

mcconsensus <- function(x,diss=FALSE,algorithms=list('agnes'),alparams=list(),clmin=2,clmax=10,prop=0.8,reps=50,numofcores=1){

  #CHECK INPUTS
  #first check if the input data is an ExpressionSet if so extract the data
  if(class(x)=='ExpressionSet'){x <- expSetProcess(x);};
  #check integrity of data (what about NAs)
  if(data_check(x)!=TRUE){stop('the provided data fails the data integrity check')}
  #check that the algorithms are passed as a list
  if(class(algorithms)!='list'){stop('You must pass the algorithm names in a list')};
  #check that all of the algorithms specified are available
  for(i in 1:length(algorithms)){
    if(existsFunction(algorithms[[i]])!=1){stop('One or more of the specified algorithms does not exist')} else {}
  }

  alweights = list()
  merge = 0

  #clmin,clmax,reps integers
  if(clmin < 2){stop('cannot have less than two clusters')}
  #prop between 0 and 1
  if(prop <= 0 | prop >1){stop('resampling proportion must be greater than 0 and less than 1')}
  #check merge


  #if we get past all this now call the re-sampling
  sample_number = dim(x)[1];

  #list to hold all of the consmatrix objects
  cmlist <- list();

  doParallel::registerDoParallel(cores=numofcores)
  #cmlistret <- apply_function(clmin:clmax, clst=cmlist);
  result <- foreach::foreach(iter=clmin:clmax) %dopar% {
    clst <- list()
    #normalise algorithm weightings
    total = sum(as.numeric(alweights));
    norm = as.numeric(alweights)/total;
    if(merge>0){
      mm <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)))
    }

    #for each algorithm
    for(a in 1:length(algorithms)){
      final_conn <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)));
      final_id <- final_conn;

      #get the algorithm to call for cm
      algo_clmem <- paste(algorithms[[a]],'_clmem',sep='');
      current_algo <- get(algo_clmem);
      #check if params have been specified
      if(length(alparams)!=0){
        current_params <- alparams[[a]];
        #check if the diss parameter is set and true, in case user forgets diss flag
        if(!is.null(current_params$diss)){
          diss=TRUE;
        }
      }
      else{
        #here someone has specified that x is a distance matrix, but not included it in the params
        if(diss==TRUE){
          current_params = list(diss=TRUE);
        }
        #otherwise x not a distance matrix params stay empty
        else{
          current_params = list();
        }
      }
      for(i in 1:reps){
        #the generic clustering algorithm calling method
        #checking if data is in fact a distance matrix (data.frame)
        #either by flag diss being set or diss being set in parameters
        if (diss==TRUE){
          #we have already put the object through the data_check
          #we need to check that the row and column names are the same in the distance data.frame
          if(unique(names(x)==row.names(x))==FALSE){stop('the row and column names of the distance matrix are not the same or the distance matrix is not square !')}
          #if OK perform the sample
          else{
            samp_row_names <- sample(row.names(x),as.integer(prop*sample_number));
            samp_x <- x[samp_row_names,samp_row_names];
            #convert to object of class 'dissimilarity'
            #samp_x <- as.dist(samp_x);
          }
        }
        else{
          samp_x <- x[sample(row.names(x),as.integer(prop*sample_number)),]; #ORIGINAL
        }
        clmem <- current_algo(samp_x,iter,params=current_params);
        conn <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(c(row.names(x)),c(row.names(x))))

        conn[] <- 0
        for(cluster in unique(clmem$cm)) {
          conn[rownames(subset(clmem,cm==cluster)),rownames(subset(clmem,cm==cluster))] <- 1
        }

        id = conn;
        id[] <- 0;
        id[row.names(clmem),row.names(clmem)] <- 1;
        final_conn = final_conn+conn;
        final_id = final_id+id;
      }
      #this is the consensus matrix
      consensus = final_conn/final_id;
      #check for NAs (this is safer than pseudo counting, replace NAs with 0 i.e. they are never drawn together and/or never connected)
      ind <- is.na(consensus)
      consensus[ind] <- 0
      #perform the full clustering as reference
      rm <- current_algo(x,iter,params=current_params)


      #now create the S4 class
      current_consmatrix <- new('consmatrix',cm=as.matrix(consensus),rm=rm,a=algorithms[[a]],k=iter);
      #add the consmatrix object to the list
      clst[paste('e',a,'_',algorithms[[a]],'_k',iter,sep='')] <- current_consmatrix
      #now add to the running matrix for merge weighting by the correct algorithm weighting if merge !=0

    }



    return (clst)
  }

  return(unlist(result));
}



#' Consboost function which using consensus clustering method
#' @param object the FCS exprs data set after transformed
#' @param K Number of final clusters
#' @param sample_n sampling the data for consensus clustering
#' @param n_core CPU thread use
#' @import ada
#' @export

cytosee_consboost <- function(object,K,sample_n = 3000,n_core=3){
  set.seed(443)
  data=object@fcs.data[object@event.use,][,object@channel.use]
  colnames(data)<-paste0("marker",1:ncol(data))
  total_len = length(data[,1])
  if(total_len<=sample_n){
    trainset <-data.class()[1:total_len,]
  }

  else{
    mm <- kmeans(data,100,iter.max=1000000)
    label <- mm$cluster
    av_num = sample_n / total_len
    label_level = levels(as.factor(label))
    return_list=vector()
    for(i in label_level){
      ll = which(label==i)
      return_list = c(return_list,sample(ll,round(length(ll)*av_num)))
    }
    trainset <- data[base::sort(return_list),]
    message("calculating trainset")
    cons_data <- mcconsensus(as.data.frame(trainset),diss=FALSE,algorithms=list('hclust'),alparams=list(),clmin=K-1,clmax=K+1,prop=0.8,reps=100,numofcores=n_core)
    trainClass <- as.vector(paste0("X",cons_data[[2]]@rm$cm))
    message("Adaboost model")
    AdaModels <- oneVsAll(trainset, trainClass, ada)
    ada.pred.unknown   <- predict(AdaModels, as.data.frame(data),  type='probs')
    Ada.pred.unknown   <- data.frame(lapply(ada.pred.unknown, function(x) x[,2]))
    ADA.pred.unknown   <- classify(Ada.pred.unknown)
    ada.class <- as.integer(gsub("X","",ADA.pred.unknown$Class))
    object@ClusterID <- as.data.frame(ada.class)
    return(object)
  }
}
