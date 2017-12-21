#' build Divisive marker tree
#'
#' description
#' This is a function for building divisive marker tree. It can be used to divide the the cell populations into their own clusters by surface markers.
#' @param data this data should be a dataframe which contain the all elements from the expression data of .fcs file and it's first column should be the cluster num of this event.
#' @example
#' "to be continue!"
#' @return DMT TREE list which is a tree list contain elements like "label,size,id,previd,clus"
#' @export


#get sd value
DMT=function(data){
  f=data
  cl_average_vec=data.frame()
  cl=unique(f[,1])
  Colnames=colnames(f[,2:length(f)])
  data=f[,Colnames]
  len_col=length(Colnames)
  for(i in c(1:length(cl))){
    temp=data[which(f[,1]==cl[i]),]
    for(j in c(1:length(Colnames))){
      cl_average_vec[i,j]=mean(as.numeric(temp[,j]))
    }
  }
  rownames(cl_average_vec)=cl
  colnames(cl_average_vec)=Colnames
  cl_average_rownames=row.names(cl_average_vec)
  sdVec=rep(0,len_col)
  for(i in c(1:length(cl))){
    temp=data[which(f[,1]==cl[i]),]
    len_temp=length(temp[,1])
    for(j in c(1:len_col)){
      sdVec[j]=sdVec[j]+var(temp[,j])*(len_temp-1)
    }
  }

  for(i in c(1:len_col)){
    sdVec[i]=max(0.25,sqrt(sdVec[i]/length(f[,1])))
  }

  #a method for dividing a group of cells into positive and negative parts by one marker
  split=function(node){
    if(length(node$clus)==2){
      cl1=node$clus[1]
      cl2=node$clus[2]
      vec1=cl_average_vec[which(row.names(cl_average_vec)==cl1),]
      vec2=cl_average_vec[which(row.names(cl_average_vec)==cl2),]
      diff=vec1-vec2
      maxCol=0
      maxdiff=0
      breakpoint=""
      diff=diff/sdVec
      for(i in c(1:length(diff))){
        diff[i]=abs(diff[i])
        diff[i]=diff[i]/sdVec[i]
        if(diff[i]>maxdiff){
          maxDiff = diff[i];
          maxCol = i;
          breakpoint = as.character(round((vec1[i]+vec2[i])/2,2))
        }
      }
      if(diff[maxCol]>0){
        label1=paste0(Colnames[maxCol],">",breakpoint)
        label2=paste0(Colnames[maxCol],"<",breakpoint)
      }
      else{
        label1=paste0(Colnames[maxCol],"<",breakpoint)
        label2=paste0(Colnames[maxCol],">",breakpoint)
      }
      node1=list(label=label1,
                 size=length(f[which(f[,1]==cl1),1]),
                 id=id+1,
                 previd=node$id,
                 clus=cl1)
      node2=list(label=label2,
                 size=length(f[which(f[,1]==cl2),1]),
                 id=id+2,
                 previd=node$id,
                 clus=cl2)
      id<<-id+2
      return(list(node1,node2))
    }
    bestCol=-1
    minSumSq=1.7976931348623157E308
    bestDivVal=-1
    bestSepRange=0
    len_clu=length(node$clus)
    for(i in c(1:len_col)){
      for(j in c(1:len_clu)){
        divVal=cl_average_vec[which(row.names(cl_average_vec)==node$clus[j]),i]
        div_group=divideClusters(node=node,col=i,divVal = divVal)
        if(length(div_group)>0){
          sumSq=getSumAngDist(div_group[[1]])+getSumAngDist(div_group[[2]])
          if(sumSq<minSumSq){
            bestDivVal=divVal
            bestCol=i
            minSumSq=sumSq
            minSepRange=1.7976931348623157E308
            for(k in c(1:length(div_group[[1]]))){
              for(l in c(1:length(div_group[[2]]))){
                sepR=abs(cl_average_vec[which(row.names(cl_average_vec)==div_group[[1]][k]),i]-cl_average_vec[which(row.names(cl_average_vec)==div_group[[2]][l]),i])/sdVec[i]
                if(sepR<minSepRange){
                  minSepRange=sepR
                }
              }
            }
            bestSepRange = minSepRange
          }
          else if(sumSq/minSumSq<1.0001){
            minSepRange=1.7976931348623157E308
            for(k in c(1:length(div_group[[1]]))){
              for(l in c(1:length(div_group[[2]]))){
                sepR=abs(cl_average_vec[which(row.names(cl_average_vec)==div_group[[1]][k]),i]-cl_average_vec[which(row.names(cl_average_vec)==div_group[[2]][l]),i])/sdVec[i]
                if(sepR<minSepRange){
                  minSepRange=sepR
                }
              }
            }
            if (minSepRange > bestSepRange) {
              bestDivVal = divVal;
              bestCol = i;
              minSumSq = sumSq;
              bestSepRange = minSepRange;
            }
          }
        }
      }
    }
    sep=divideClusters(node,bestCol,bestDivVal)
    avg1=mean(f[which(f[,1] %in% sep[[1]]),bestCol+1])
    avg2=mean(f[which(f[,1] %in% sep[[2]]),bestCol+1])
    val1=c()
    val2=c()

    for(i in c(1:length(sep[[1]]))){
      val1[i] = cl_average_vec[which(row.names(cl_average_vec)==sep[[1]][i]),bestCol]
    }
    for(i in c(1:length(sep[[2]]))){
      val2[i] = cl_average_vec[which(row.names(cl_average_vec)==sep[[2]][i]),bestCol]
    }

    val1_max=max(val1)
    val1_min=min(val1)
    val2_max=max(val2)
    val2_min=min(val2)

    if(avg1>=avg2){
      breakpoint=as.character(round((val1_max+val2_min)/2,2))
      label1=paste0(Colnames[bestCol],">",breakpoint)
      label2=paste0(Colnames[bestCol],"<",breakpoint)
    }
    else{
      breakpoint=as.character(round((val2_max+val1_min)/2,2))
      label1=paste0(Colnames[bestCol],"<",breakpoint)
      label2=paste0(Colnames[bestCol],">",breakpoint)
    }
    node1=list(label=label1,
               size=length(f[which(f[,1] %in% sep[[1]]),1]),
               id=id+1,
               previd=node$id,
               clus=sep[[1]])
    node2=list(label=label2,
               size=length(f[which(f[,1]%in% sep[[2]]),1]),
               id=id+2,
               previd=node$id,
               clus=sep[[2]])
    id<<-id+2
    return(list(node1,node2))
  }

  #a function to divede vectors into two parts
  divideClusters=function(node,col,divVal){
    g1=c()
    g2=c()
    for(i in c(1:length(node$clus))){
      if(cl_average_vec[which(row.names(cl_average_vec)==node$clus[i]),col]<divVal){
        g1=c(g1,node$clus[i])
      }
      else{
        g2=c(g2,node$clus[i])
      }
    }
    if(length(g1)==0 ||length(g2)==0 ){
      return()
    }
    return(list(g1,g2))
  }

  #function for getting sum Distance
  getSumAngDist=function(cl){
    cl_aver_vec_sd=c()
    cl_row=which(cl_average_rownames %in% cl)
    for(i in c(1:length(Colnames))){
      cl_aver_vec_sd[i]=mean(cl_average_vec[cl_row,i])/sdVec[i]
    }
    sumDist=0
    vec=cl_average_vec[cl_row,]
    vec=as.matrix(vec)
    for(i in c(1:length(cl))){
      for(j in c(1:length(sdVec))){
        vec[i,j]=vec[i,j]/sdVec[j]
      }
    }
    for(c in c(1:length(cl))){
      sumDist=sumDist+acos(min(1,sum(vec[c,]*cl_aver_vec_sd)/(sqrt(sum(vec[c,]*vec[c,]))*sqrt(sum(cl_aver_vec_sd*cl_aver_vec_sd)))))
    }
    return(sumDist)
  }


  #Build DMT_TREE

  root=list(label="Root",
            size=length(f[,1]),
            id=0,
            previd=NULL,
            clus=as.character(cl)
  )
  DMT1=list()
  DMT1[[1]]=root
  last_step=list()
  last_step[[1]]=root
  id=0
  while(length(last_step)>0){
    next_step=list()
    for(node in last_step){
      children=split(node)
      for(ch in children){
        DMT1=c(DMT1,list(ch))
        if(length(ch$clus)>1){
          next_step=c(next_step,list(ch))
        }
      }
    }
    last_step=list()
    last_step=next_step
  }
  return(DMT1)
}



##### get avg vec #####

### inner method not for users ###
getAvg=function(fcsData,ClusterID){
  data=cbind(ClusterID,fcsData)
  f=data
  cl_average_vec=data.frame()
  cl=unique(f[,1])
  Colnames=colnames(f[,2:length(f)])
  data=f[,Colnames]
  len_col=length(Colnames)
  for(i in c(1:length(cl))){
    temp=data[which(f[,1]==cl[i]),]
    for(j in c(1:length(Colnames))){
      cl_average_vec[i,j]=mean(temp[,j])
    }
  }
  rownames(cl_average_vec)=cl
  colnames(cl_average_vec)=Colnames
  return(cl_average_vec)
}

