### Prepared by Zaid Abdo/Metagenomics Class-2017
######## CRAN packages
#### The next three lines need to all be un hashed to install packages from CRAN
#### that you need to load using the library() function. This should be done
#### only once the first time you install R or when you upgrade it to a new 
#### vergion

#install.packages(c("pheatmap","vegan","ggplot2","cluster","NbClust"
#                 ,"gplots","devtools","RColorBrewer","corrgram",
#                 "fdrtool","vegan3d","rgl"))

#install.packages(c("pheatmap","vegan","ggplot2","gplots","corrgram","vegan3d","rgl"))


#### Use library() to load the packages you installed from CRAN
## pheatmap: creates heatmaps with clustering package and is used in the 
## functions pheatmap.1.ftn(), pheatmap.2.ftn() and pheatmap.3.ftn
library(pheatmap)
## vegan: is an evological analysis package 
library(vegan)
## ggplot2: a plotting package used for the bar plots functions: bar.taxa.sample.ftn() and bar.sample.trt.ftn()
library(ggplot2)
## plotting package used for heatmaps (function heatmap.2 in Analysis-Normalization-Ordination.R)
library(gplots)
## grid: A graphic package required to fine tune some of the other packages
library(grid)
## vegan3d and rgl: required packages for 3D ordination plots using vegan 
## in Analysis-Normalizatiopn-Ordination.R (see description there)
library(vegan3d)
library(rgl)

######## Bioconductor
#### Use source to install bioclite.R which will allow you to load packages from Bioconductor
#source("https://bioconductor.org/biocLite.R")
#### Run biocLite() to install dependencies required to install other packages. This should be done only once when you inastell
#### or when you update R.
#biocLite()
#### Run biocLite to install the packages: phyloseq, metagenomSeq and DESeq2 that we use for data management
#### and data analysis
#biocLite('phyloseq')
#biocLite('metagenomeSeq')
#biocLite('DESeq2')

#### use library() to load these packages
## phyloseq: we use this in all .R programs for data management and other analyses
library(phyloseq)
## metagenomeSeq: we use this to perform the Cumulative Sum Scaling normalization in Analysis-Normalization-Ordination.R and
## for data analysis in Analysis-Statistics.R
library(metagenomeSeq)
## DESeq2 we use this in Analysis-Statistics.R for another way to perform data nalaysis using an approach from RNA-Seq
library(DESeq2)