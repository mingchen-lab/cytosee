## CytoSEE: A Web-based toolkit for automatic computation and evaluation of cytometry data


### Workflow of CytoSEE
![CytoSEE](http://bis.zju.edu.cn/picture/workflow_new_cytosee.png) 





### Installation
Currently, A docker version is available at dockerhub [cytosee](https://hub.docker.com/r/mchenlab/cytosee), you can also edit and build your own image with [cytosee Dockerfile]()

#### Linux 

```R
install.package("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite(c('FlowSOM','flowMeans','SamSPECTRAL'),ask=FALSE)
devtools::install_github('madlogos/recharts')
devtools::install_github('mingchen-lab/cytosee')
```

#### Windows 

Rtool34 or higher version is needed depended on the R version. Download address: https://www.r-project.org/
```R
install.package("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite(c('FlowSOM','flowMeans','SamSPECTRAL'),ask=FALSE)
devtools::install_github('madlogos/recharts')
devtools::install_github('mingchen-lab/cytosee')
```


### Usage 
```R
# CytoSEE start with one command in your IDE(R/Rstudio)
library("cytosee")
cytosee_gui()
```


### NEWS

We have integrate [FIt-SNE](https://github.com/KlugerLab/FIt-SNE) into CytoSEE to deal with dimention reduction of large cell counts. To enable this function, FIt-SNE should be avaliable in your environment.


