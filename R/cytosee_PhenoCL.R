#' Automate Query phenotypes from Cell ontology
#'
#' This function will automatically produce Phenotypes for each clusters and query the most possible result for them.
#' @name cytosee_PhenoCL
#' @aliases cytosee_PhenoCL
#' @param data A matrix contain fcs data .
#' @param ClusterID A data.frame contain the Cluster ID of each events in data.
#' @param MarkerName character vectors for revising all the name of markers. Default is a empty matrix.
#' @param MaxLabelnum A numbmic vector determine how many top labels will be display. Default is 10.
#' @return A list of all most possible results of clusters.
#' @export
#' @examples
#' # load the library
#' library(cytosee)
#' data(GvHD001)
#' data=GvHD001$data[,c(3:6)]
#' ClusterID=GvHD001$label
#' MarkerName=c("CD4","CD8b","CD8","CD3")
#' result=cytosee_PhenoCL(data,ClusterID,MarkerName)
#'
cytosee_PhenoCL<-function(data,ClusterID,MarkerName=c(),MaxLabelnum=10){
  if(!is.null(MarkerName)){
    colnames(data)=MarkerName
  }
  if(length(grep(pattern = "[, .(*]",colnames(data)))>0){
    stop("Illegal symbol was detected such as '(' in your marker names")
  }
  data=cbind(ClusterID,data)
  message("Run cytosee_DMT...")
  dmtr=cytosee_DMT(data)
  markers=c()
  message("Data transformation...")

  for(i in c(1:length(dmtr)-1)){
    re=unlist(strsplit(dmtr[[i+1]]$label,"[<||>]"))
    markers[i]=re[1] #in order to delet the root
  }

  # simplify marker stucture
  markers=unique(markers)


  select_marker=list()
  Cluster_marker=c()
  Cluster_m2c=c()
  cnt=0
  for(i in c(1:length(dmtr))){
    if(length(dmtr[[i]]$clus)==1){
      tmp=dmtr[[i]]
      str=""
      while(tmp$id!=0){
        str=paste0(tmp$label,",",str)
        tmp=dmtr[[tmp$previd+1]]
      }
      cnt=cnt+1
      Cluster_marker[cnt]=str
      Cluster_m2c[cnt]=dmtr[[i]]$clus
    }
  }
  # make the phenotype eazier for transform into style of "+","-" style.
  re=gsub(pattern = " ",replacement = "",Cluster_marker)
  re=gsub(pattern = "<.*?,",replacement = "<",re)
  re=gsub(pattern = ">.*?,",replacement = ">",re)

  for(s in c(1:length(markers))){
    word=markers[s]
    word_p=paste0(word,">")
    word_m=paste0(word,"<")
    for(i in c(1:length(re))){
      store=re[i]
      # first and second symbol will determine this marker is "+" or "-"
      first=""
      second=""
      cnt=0
      while(length(grep(pattern = word_m,x = store,fixed = TRUE))!=0 || length(grep(pattern =word_p,x = store,fixed = TRUE))!=0 ){
        if(length(grep(word_p,store))!=0){
          cnt=cnt+1
          store=sub(pattern=word_p,replacement="",store,perl = FALSE)
          if(cnt==1){
            first="plus"
          }
          else if(cnt==2){
            second="plus"
          }
        }
        if(length(grep(word_m,store))!=0){
          cnt=cnt+1
          store=sub(pattern = word_m,replacement="",store,perl = FALSE)
          if(cnt==1){
            first="minus"
          }
          else if(cnt==2){
            second="minus"
          }
        }
      }
      if(cnt==0){
        next()
      }
      if(first=="plus"&&second=="plus"){
        re[i]=paste0(store,word,"++")
      }
      else if(first=="plus"){
        re[i]=paste0(store,word,"+")
      }
      else if(first=="minus"&&second=="minus"){
        re[i]=paste0(store,word,"--")
      }
      else if(first=="minus"){
        re[i]=paste0(store,word,"-")
      }
    }
  }

  # query the result from local Cell ontology
  message("Querying the label......")
  CL_label=list()
  for(i in 1:length(re)){
    result=cytosee_LocCL(MarkerList = re[i],MaxHitsPht=MaxLabelnum)
    ClusterID=Cluster_m2c[i]
    CL_label[[i]]=list("ClusterID"=ClusterID,"result"=result,"PhenoType"=re[i])
  }
  message("Done!")
  return(CL_label)
}





#' @title build Divisive marker tree(DMT)
#'
#' description
#' This is a function for building divisive marker tree. It can be used to divide the the cell populations into their own clusters by surface markers.
#' @param data this data should be a dataframe which contain the all elements from the expression data of .fcs file and it's first column should be the cluster num of this event.
#' @return DMT TREE list which is a tree list contain elements like "label,size,id,previd,clus"
#'
cytosee_DMT=function(data){
  if(!is.data.frame(data)){
    stop("data must be a dataframe!")
  }
  f<-data
  cl_average_vec=data.frame()
  cl=unique(f[,1])
  Colnames=colnames(f[,2:length(f[1,])])
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
    if(len_temp==1){
      warning("There is a cluster whose cell number is '1',we will regrad 'sd' value as zero")
    }
    for(j in c(1:len_col)){
      if(len_temp==1){
        sdVec[j]=0
      }
      else{
        sdVec[j]=sdVec[j]+var(temp[,j])*(len_temp-1)
      }
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


#' @title It will query the possible result from local ontology datatable
#'
#' @description
#'
#' @param MarkerList a string word such as "CD19-CD3+CD49b++Sca1-".
#' @param MaxHitsPht a numerical value, the number of max possible value.
#' @return if plot is TRUE ,The plot list will be returned, else a list of top 5 results and their score will be returned
#' @export
#' @examples
#'
#' MarkerList="CD8+CD4+CCR7-"
#' cytosee_LocCL(MarkerList)
#'
cytosee_LocCL <- function (MarkerList, MaxHitsPht=10){
  #surface marker divided
  markerlist=unlist(strsplit(split = "\\++|-+",MarkerList))
  expresslist=unlist(strsplit(split = "\\w+",MarkerList))
  Markertype=list(bright="",positive="",low="",negative="")
  for( i in 1:length(markerlist)){
    if(expresslist[i+1]=="++"){
      Markertype["bright"]=paste0(markerlist[i],",",Markertype$bright)
    }
    else if(expresslist[i+1]=="+"){
      Markertype["positive"]=paste0(markerlist[i],",",Markertype$positive)
    }
    else if(expresslist[i+1]=="--"){
      Markertype["low"]=paste0(markerlist[i],",",Markertype$low)
    }
    else if(expresslist[i+1]=="-"){
      Markertype["negative"]=paste0(markerlist[i],",",Markertype$negative)
    }
  }



  # marker in alias should make sense for judge whether this label is fit for our query.
  Alias_score_label=c()

  #qurey the result
  res=list()
  for(i in 1:length(Markertype)){
    if(Markertype[i]!=""){
      marker=unlist(strsplit(as.character(Markertype[i]),split=","))
      for(j in marker){
        if(i==1){
          #change the string for reg
          markername=paste0(j,"++")
          j=paste0(j,"[a-z]?")
          j=paste0(j,"-bright","|",j,"-positive","|",j,"[+]")
          name=CL_lib[which(grepl(j,CL_lib$label)),2]
          alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
          comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
          comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
          uni=Reduce(union,list(name,alias,comment,comment2))
          res[markername]=data.frame(uni)
          Alias_score_label=c(Alias_score_label,as.character(alias),as.character(name))
        }
        else if(i==2){
          markername=paste0(j,"+")
          j=paste0(j,"[a-z]?")
          j=paste0(j,"-positive","|",j,"[+]")
          name=CL_lib[which(grepl(j,CL_lib$label)),2]
          alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
          comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
          comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
          uni=Reduce(union,list(name,alias,comment,comment2))
          res[markername]=data.frame(uni)
          Alias_score_label=c(Alias_score_label,as.character(alias),as.character(name))
        }
        else if(i==3){
          markername=paste0(j,"--")
          j=paste0(j,"[a-z]?")
          j=paste0(j,"-low","|",j,"-dim","|",j,"-negative","|",j,"- negative")
          name=CL_lib[which(grepl(j,CL_lib$label)),2]
          alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
          comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
          comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
          uni=Reduce(union,list(name,alias,comment,comment2))
          res[markername]=data.frame(uni)
          Alias_score_label=c(Alias_score_label,as.character(alias),as.character(name))
        }
        else if(i==4){
          markername=paste0(j,"-")
          j=paste0(j,"[a-z]?")
          j=paste0(j,"-negative","|",j,"- negative")
          name=CL_lib[which(grepl(j,CL_lib$label)),2]
          alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
          comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
          comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
          uni=Reduce(union,list(name,alias,comment,comment2))
          res[markername]=data.frame(uni)
          Alias_score_label=c(Alias_score_label,as.character(alias),as.character(name))
        }
      }
    }
  }

  #choose the result which is most possible cell type for the query
  all_label=c()

  # get all result and proviede a score for them
  marker_score=c()
  marker_alias_score=c()
  level_score=c()
  child_score=c()
  final_score=c()

  # calculate marker score
  for(i in length(res):1){
    len_queue=length(t((combn(length(res),i)))[,1])
    for(j in 1:len_queue){
      result=Reduce(intersect,res[t(combn(length(res),i))[j,]])
      for(k in result){
        if(length(k)!=0 && !(k %in% all_label)){
          all_label=c(all_label,c(k))
          marker_score=c(marker_score,c(i))
        }
      }
    }
  }
  # calculate level score
  treelist=list()
  for( i in 1:length(all_label)){
    treelist[i]=all_label[i]
    child=all_label[i]
    rank=1
    while(!is.na(CL_lib[which(array(CL_lib[,2])==child),4])){
      parent=CL_lib[which(array(CL_lib[,2])==child),4]
      treelist[i]=paste0(treelist[i],">",parent)
      child=parent
      rank=rank+1
    }
    level_score=c(level_score,1/rank)
  }
  # calculate child score and marker_alias_score
  for(i in 1:length(all_label)){
    rank=length(grep(all_label[i],treelist,fixed = TRUE))
    child_score=c(child_score,rank)
    marker_alias_score=c(marker_alias_score,length(Alias_score_label[which(Alias_score_label==all_label[i])])+1)
  }
  # calculate final score
  for(i in 1:length(all_label)){
    rank=marker_score[i]/max(marker_score)*0.6+level_score[i]/max(level_score)*0.1+child_score[i]/max(child_score)*0.1+marker_alias_score[i]/max(marker_alias_score)*0.2
    final_score=c(final_score,rank)
  }

  #get top label
  tabel_scorelist=cbind(all_label,Score=round(final_score,3))
  tabel_scorelist=as.data.frame(tabel_scorelist[order(tabel_scorelist[,2],decreasing = TRUE),])
  top_label=head(x = tabel_scorelist,MaxHitsPht)$all_label
  top_score=head(x = tabel_scorelist,MaxHitsPht)$Score
  label_score=top_score

  if(is.null(top_label)){
    CL_result=list("QueryMarker"=MarkerList,"Labels"=top_label,"Score"=label_score,"Nodes"=NULL,"Links"=NULL)
    warning("There are some cell phenotypes can not get the query result!")
    return(CL_result)
  }

  #parent tree 1. query the parent node from CL_lib
  treelist=list()
  for( i in 1:length(top_label)){
    treelist[i]=as.character(top_label[i])
    child=top_label[i]
    rank=1
    while(!is.na(CL_lib[which(array(CL_lib[,2])==child),4])){
      parent=CL_lib[which(array(CL_lib[,2])==child),4]
      treelist[i]=paste0(treelist[i],">",parent)
      child=parent
    }
  }

  #parent tree 2.build Links
  len_cells=0
  for(i in 1:length(treelist)){
    temp=unlist(strsplit(treelist[[i]],">"))
    len_cells=len_cells+length(temp)
  }
  from=c()
  to=c()
  for(i in 1:length(treelist)){
    tmplist=list()
    tmplist[i]=matrix(ncol = length(treelist[i]))
    tmplist[i]=strsplit(unlist(treelist[i]),">")
    for(j in 1:length(tmplist[[i]])){
      if(j==length(tmplist[[i]])||tmplist[[i]][j]=="native cell"){
        next()
      }
      else{
        from=c(tmplist[[i]][j],from)
        to=c(tmplist[[i]][j+1],to)
      }
    }
  }
  Links=data.frame(from,to)
  Links=as.data.frame(unique(Links))

  ### Nodes paras for drawing Tree###
  id=union(unique(from),unique(to))
  color=c()
  title=c()
  for(i in id){
    if(i %in% top_label){
      color=c(color,"#BFDBB8")
      t="Markers:"
      for(j in names(res)){
        re=eval(parse(text=sprintf("res$'%s'",j)))
        if(i %in% re){
          t=paste0(t,j)
        }
      }
      title=c(title,t)
    }
    else{
      color=c(color,"#DBCBB8")
      title=c(title,"None")
    }
  }

  Nodes=data.frame(id=id,label=id,color=color,title=title)
  CL_result=list("QueryMarker"=MarkerList,"Labels"=top_label,"Score"=label_score,"Nodes"=Nodes,"Links"=Links)
  return(CL_result)
}






