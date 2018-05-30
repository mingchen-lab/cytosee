## CytoSEE: R based toolkit for automatic computation and evaluation of cytometry data


### Workflow of CytoSEE
![CytoSEE](http://bis.zju.edu.cn/picture/workflow_new_cytosee.png) 





### Installation

#### Linux Ubuntu 

some package will rely on these tools:
 
```Shell
apt-get install -y wget r-base  supervisor \
libcurl4-openssl-dev
libssl-dev 
libcgal-dev 
libglu1-mesa-dev
libglu1-mesa-dev 
libcairo2-dev
libxt-dev
gdebi-core 
pandoc
pandoc-citeproc
libxml2-dev 
```

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


