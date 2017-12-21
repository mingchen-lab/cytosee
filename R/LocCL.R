#' It will query the possible result from local ontology file
#' description
#' @param MarkerList a string word such as "CD19-CD3+CD49b++Sca1-".
#' @param MaxHitsPht a numerical value, the number of max possible value.
#' @param plot a boolen value , It's not suitable for user to choose, default is FALSE.
#' @return if plot is TRUE ,The plot list will be returned, else a list of top 5 results and their score will be returned
#' @export LocCL
#' @example
#' LocCL()



LocCL <- function ( MarkerList = "CD19-CD3-CD49b+Sca1-CD150+CD16-", MaxHitsPht=5,plot="FALSE"){
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
            j=paste0(j,"-bright","|",j,"-positive","|",j,"\\+")
            alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
            comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
            comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
            uni=Reduce(union,list(alias,comment,comment2))
            res[markername]=data.frame(uni)
          }
          else if(i==2){
            markername=paste0(j,"+")
            j=paste0(j,"[a-z]?")
            j=paste0(j,"-positive","|",j,"\\+")
            alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
            comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
            comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
            uni=Reduce(union,list(alias,comment,comment2))
            res[markername]=data.frame(uni)
          }
          else if(i==3){
            markername=paste0(j,"--")
            j=paste0(j,"[a-z]?")
            j=paste0(j,"-low","|",j,"-dim","|",j,"-negative","|",j,"- negative")
            alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
            comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
            comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
            uni=Reduce(union,list(alias,comment,comment2))
            res[markername]=data.frame(uni)
          }
          else if(i==4){
            markername=paste0(j,"-")
            j=paste0(j,"[a-z]?")
            j=paste0(j,"-negative","|",j,"- negative")
            alias=CL_lib[which(grepl(j,CL_lib$alias)),2]
            comment=CL_lib[which(grepl(j,CL_lib$comment)),2]
            comment2=CL_lib[which(grepl(j,CL_lib$comment2)),2]
            uni=Reduce(union,list(alias,comment,comment2))
            res[markername]=data.frame(uni)
          }
        }
      }
    }

    #choose the result which is most possible cell type for the query
    hitsPht=1
    top_label=c()

    #score the result
    label_score=c()
    for(i in length(res):1){
      len_queue=length(t((combn(length(res),i)))[,1])
      for(j in 1:len_queue){
        result=Reduce(intersect,res[t(combn(length(res),i))[j,]])
        for(k in result){
          if(length(k)!=0 && !(k %in% top_label)){
            top_label=c(top_label,c(k))
            label_score=c(label_score,c(i))
          }
          if(length(top_label)>=MaxHitsPht){
            break
          }
        }
        if(length(top_label)>=MaxHitsPht){
          break
        }
      }
      if(length(top_label)>=MaxHitsPht){
        break
      }
    }

    #parent tree 1. query the parent node from CL_lib
    treelist=list()
    for( i in 1:length(top_label)){
      treelist[i]=top_label[i]
      child=top_label[i]
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

    ### nodes paras ###
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
    CL_result=list("Labels"=top_label,"Score"=label_score,"Nodes"=Nodes,"Links"=Links)
    return(CL_result)
}


