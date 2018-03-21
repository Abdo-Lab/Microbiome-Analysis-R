######## Functions prepared by Zaid Abdo - 12/2017 - Microbial Genomics/Metagenomics Data Analysis
### Species bar plot per sample using most abundant in sample
## This function requires an otu table (on any taxonomic level) in the form of a data frame
## It only plots the most abundant based on the cutoff given (the cutoff defaults to 0 and hence plots everything)
## x label defaults to "Sample" but can be changed, and y label defaults to "OTU" but can be changed
## The function plots proportions and it is best to use the non-normalized data for this. It doesn't split the plot
## per treatment level; but combines all samples together.
bar.taxa.sample.ftn = function(d.df,cutoff=0,xlab="Sample",llab="OTU"){
  d.rs = apply(d.df,1,sum)
  ln = length(d.rs)
  p.df = d.df/d.rs
  if(cutoff > 0){
    cut.b = apply(p.df,2,function(x)sum(x>cutoff))
    p.df = p.df[,cut.b > 0]
    p.rs = apply(p.df,1,sum)
    p.df = p.df/p.rs
  }
  nm = colnames(p.df)
  r.nm = row.names(p.df)
  sp = rep(nm,ln)
  sample = w = c()
  for(i in 1:ln){
    sample = c(sample,rep(r.nm[i],length(nm)))
    w = c(w,t(p.df[i,]))
  }
  hist.df = data.frame(sample,sp,w)
  print(qplot(data=hist.df,x=sample,fill=sp,geom="bar",weight=w)
        +geom_bar(color="white",linetype=1)
        +labs(x=xlab,y="Proportion")+scale_fill_discrete(name=llab)
        +theme(axis.text.x=element_text(angle=0))
        +scale_color_manual(values = palette(rainbow(19)))
        +guides(fill=guide_legend(ncol=1))
        +theme_minimal())

}

### Species bar plot per sample using most abundant in sample split by treatment level
## This function requires an otu table (on any taxonomic level) in the form of a data frame
## It only plots the most abundant based on the cutoff given (the cutoff defaults to 0.1)
## title defaults to "Experiment", x label defaults to "Sample" but can be changed, and 
## y label defaults to "OTU" but can be changed. The function plots proportions and it is
## best to use the non-normalized data for this. 
## It does splits the plot per treatment level as provided by Trt, which is required. Trt is a data frame 
## that includes either 2 columns (the first is a treatment level and the second is the sample names) or
## 3 columns (first two columns include treatment levels for two factors and last includes sample names)
bar.sample.trt.ftn = function(data,Trt, cutoff=0.01, title="Experiment",xlab="Sample",llab="OTU"){
  ln = length(data[,1])
  p.df = data/apply(data,1,sum)
  nm = colnames(p.df)
  if(cutoff > 0){
    cut.b = apply(p.df,2,function(x)sum(x>cutoff))
    p.df = p.df[,cut.b > 0]
    p.rs = apply(p.df,1,sum)
    p.df = p.df/p.rs
  }
  nm.trt = colnames(Trt)
  hist.df = sp = w = c()
  for(j in 1:length(p.df[1,])){
    sp = c(sp,rep(nm[j],ln))
    w = c(w,p.df[,j])
    hist.df = rbind(hist.df,Trt)
  }
  hist.df = data.frame(hist.df,sp,w)
  names(hist.df) = c(nm.trt,"sp","w")
  if(length(hist.df[1,]) == 4){
    ### one treatment + samples
    hist.df = hist.df[order(hist.df[,1],hist.df[,2]),]
    print(ggplot(data=hist.df,aes(x=hist.df[,2],fill = hist.df[,3]))
          +geom_bar(aes(weight=hist.df[,4]),color="white",linetype=1)
          +facet_wrap(~hist.df[,1],scales = "free")
          +labs(title=title,colour=llab,x=xlab,y="Proportion")
          +scale_fill_discrete(name = llab)          
          +theme_minimal()
          +guides(fill=guide_legend(ncol=1))
          +theme(axis.text.x=element_text(angle=90)))
  }else{
    ### two treatments + samples
    hist.df = hist.df[order(hist.df[,1],hist.df[,2],hist.df[,3]),]
    print(ggplot(data=hist.df,aes(x=hist.df[,3],fill=hist.df[,4]))
          +geom_bar(aes(weight=hist.df[,5]),color="white",linetype=1)
          +facet_wrap(hist.df[,1]~hist.df[,2],scales = "free")
          +labs(title=title,colour=llab,x=xlab,y="Proportion")
          +scale_fill_discrete(name = llab)          
          +theme_minimal()
          +guides(fill=guide_legend(ncol=1))
          +theme(axis.text.x=element_text(angle=90)))
  }
}

######### Clustering and heatmaps using pheatmap
### This function plots heatmaps and clusters the different samples identifying the different treatment levels
### per sample. It requires an otu table (on any taxonomic level) in a data frame format, and the experimental
### design (meta.df data frame) with sample names as last column and either 1 factor or two factors (meta.df 
### with either 2 or three columns). Number of clusters is set to three by defauls but can be changed, in will 
### cluster by column if col.clust = T using euclidian distance but col.clust can equal to a clustering using 
### the function hclust (see Analysis-Normalization-Ordination-Clustering.R file). Cell height and width default
### to 10 but can be changed and font size defaults to 10 but can be changed.
pheatmap.3.ftn = function(d.df,trt.df, cutoff=3,col.clust = T, cellh=10,cellw=10,fonts=10){
  if(length(trt.df)==2){
    Trt = trt.df[,1,drop=FALSE]
    print(pheatmap(as.matrix(d.df),
                   cluster_rows = F,
                   cluster_cols = col.clust,
                   cutree_rows = cutoff,
                   cellwidth = cellw,
                   cellheight = cellh,
                   annotation_col = Trt,
                   fontsize = fonts,
                   height = 7,
                   width = 6.5))
  }else{
    Trt = trt.df[,1:2]
    print(pheatmap(as.matrix(d.df),
                   cluster_rows = F,
                   cluster_cols = col.clust,
                   cutree_rows = cutoff,
                   cellwidth = cellw,
                   cellheight = cellh,
                   annotation_col = Trt,
                   fontsize = fonts,
                   height = 7,
                   width = 6.5))
  }
}

###### Non metric multidimensional scalling ordination plot
#### This function creates an ordination plot using an otu table provided by the user and a vector of treatment levels 
#### that is formated as factor. The otu table can be provided on any taxonomic level. The title defaults to "Experiment"
#### but can be changed, the dimensions of the plot are set by x and y (default between -1 and 1 but should be changed)
#### there are 7 types that can be used for plotting (numbers 1-7). font size can be set by cex and line width by lwd.
#### distance used in ordination can be set by distance and defaults to "bray" (bray-curtis)
#### plot types: 1) elipsoides colored by treatment, 2) only sample names with color by treatment,
#### 3) spider (ray) plot colored by treatment 4) elipsoide and sample names, 5) spider and sample names,
#### 6) elipsoid and spider, and 7) elipsoide, spider and sample names.
ord.plot.ftn <- function(d.df,Trt,title="Experiment",x=c(-1,1),y=c(-1,1),type=1,cex=0.5,lwd=1,distance="bray"){
  l = levels(factor(Trt))
  c = seq(1,length(l))
  c1 = c()
  for(i in 1:length(l)){
    c1 = c(c1,rep(i,length(Trt[Trt==l[i]])))
  }
  s <- apply(d.df,2,sum)
  d.df <- d.df[,s>0]
  ord <- metaMDS(d.df,distance=distance,try=100,trymax=10000)
  plot(x=NULL,y=NULL,xlim=x,ylim=y,xlab="NMDS1", ylab="NMDS2")
  title(main = title)
  if(type==1){
    ordiellipse(ord, Trt, col=c,lwd=lwd,label=TRUE,cex=cex)
  }else if(type == 2){
    text(ord,display="site",cex=cex,col=c1) 
  }else if(type == 3){
    ordispider(ord, Trt, col=c, cex = cex, label = TRUE)
  }else if(type == 4){
    text(ord,display="site",cex=cex,col=c1) 
    ordiellipse(ord, Trt, col=c,lwd=lwd,label=TRUE)
  }else if(type == 5){
    text(ord,display="site",cex=cex,col=c1) 
    ordispider(ord, Trt, col=c, cex = cex, label = TRUE)
  }else if(type==6){
    ordiellipse(ord, Trt, col=c,lwd=lwd,label=TRUE)
    ordispider(ord, Trt, col=c, cex = cex, label = FALSE)
  }else{
    text(ord,display="site",cex=cex,col=c1) 
    ordiellipse(ord, Trt, col=c,lwd=lwd,label=TRUE)
    ordispider(ord, Trt, col=c, cex = cex, label = FALSE)
  }
}

