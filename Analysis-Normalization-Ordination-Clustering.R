####### If you are running the Data processing step from here (source(file="Data_Processing.R) below)
#### then running the libraries.R (source(file="libraries.R")) is redundant but doesn't hurt
source(file="libraries.R")
#### We will use some of the functions I programmed for you (have a look in the functions.R file to see what they do)
source(file="functions.R")

######## Data Processing
#### If you are happy with your data processing step you can run it again here from the "Data_Processing.R" file without
#### having to open the file again
source(file = "Data_Processing.R")

######## We use normalized data to do ordination and other analyses
##### We normalize using the Cumulative Sum Scaling (CSS) of metagenome Seq and to do so we need to create
##### a metagenomeSeq experiment using the function newMRexperiment()
##### Before that we have to do some house keeping by changing the format of the meta.df (the experimental design)
##### and the taxa.df (the taxonomic assignment) to a formate acceptable by metagenomeSeq (we use AnnotateDataFtame() function)
trt.an = AnnotatedDataFrame(meta.df)
taxa.an = AnnotatedDataFrame(taxa.df)
##### we can now create a metagenomeSeq container using newMRexperiment using the transposed raw data (t(data.df))
##### and the newly formated experimental design (trt.an) and taxonomic data (tax.an)
d.ms = newMRexperiment(t(data.df),phenoData = trt.an,featureData = taxa.an)
##### To see what is in the new metagenomeSeq container
##### We use pData() to look at the expeimental design
pData(d.ms)
##### We use fData to look at the taxonomic assignment
fData(d.ms)
##### We use MRcount() to look at the otu table (here head() shows the first 6 lineas of that table)
head(MRcounts(d.ms))

##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percintile that we will use to 
### normalize based on and we add the result in variable p
p = cumNormStatFast(d.ms)
### We adjust our container to contain the normalization pecentile using cumNorm()
d.ms = cumNorm(d.ms,p)
### We find the normalization factors per sample using normFactors()
nf = normFactors(d.ms)
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(d.ms)
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)

###### Nonmetric Multidimesional Scaling (NMDS)
#### We want to create an NMDS ordination plot and distinguish samples that belong to the different treatment levels (A and B)
#### We extract the treatment levels (column 1 from meta.df) and put this in a new vector called Trt
Trt = meta.df$trt
#### We make sure that these treatment levels are factors
Trt = factor(Trt)
#### We create a new .pdf file
pdf(file="partial-ordiplot-otu-norm.pdf",width=10,height=10)
#### We use function ord.plot.ftn() found in functions.R to plot the NMDS which needs vegan to do so
#### based on the normalized data (d.df = norm.mt), using the treatment levels (Trt), with title "NMDS plot",
#### x axes ranges from -0.5 to 0.5 [you need to play around with this to make sure the plot fits in the panel],
#### y axes ranges from -0.5 to 0.5 [you need to play around with this to make sure the plot fits in the panel],
#### plot type is number 7 (see function in functions.R),
#### font size 1 (cex=1), line width 2 (lwd=2) and distance measure is bray-curtis (distance = "bray")
#### you can change type, cex and lwd to make the plot look nice
ord.plot.ftn(d.df=norm.mt,Trt=Trt,title="NMDS plot",x=c(-0.5,0.5),y=c(-0.5,0.5),type=3,cex=1,lwd=2,distance="bray")
#### save
dev.off()

#### For a 3D NMDS
#### we first ordinate using NMDS (using function metaMDS from vegan) the normalized data (norm.mt), using bray-curtis distance
#### (distance ="bray"), we need 3 new dimensions (k = 3), and we iteriate to optimize our fit between 100 to 10000 times
ord <- metaMDS(norm.mt,distance="bray",k=3,try=100,trymax=10000)
#### now we use vegan3d and rgl packages (and xQuartz)
#### We open a 3D panel (an empty one)
open3d()
#### We plot the data points coming from the ord (the ordination object), with x, y, and z axes in black (ax.col="black)
#### and plot the points using type text (type = "t") and tips that describe the treatment levels (text = c(rep("A",5),rep("B",4)))
ordirgl(ord,ax.col="black",type="t",text= meta.df$trt)
#### We try and plot an elipsoide per treatment level (groups=Trt) and using the standard error (kind ="se")
#### with confidence region at 0.9 (conf = 0.9), and we shade (type = "shade") with colors "red" and "blue"
#### and make the shades transparent (not solid) with alpha = 0.25
orglellipse(ord,groups=Trt, kind = "se", conf = 0.9,type="shade",col=c("red","blue"),alpha=0.25)
#### We also plot spider or ray lines centered on the centroid of the 3D plot using ordspider()
orglspider(ord,groups=Trt,col=c("red","blue"))
#### Saving the plot you like
rgl.snapshot(filename="figures/ordiplot-3d-fecal-all.png", fmt = "png", top = TRUE )


####### Creating heatmaps and clustering on the otu level
#### We first calucalte a distance matrix that we use for the clustering using vegan function vegdist() on
#### the normalized data (norm.mt) and using distance measure bray-curtis (method="bray")
data.dist <- vegdist(norm.mt, method = "bray")
#### then we use the distance matrix to cluster using hclust() function and average linkage (UPGMA)
row.clus <- hclust(data.dist, "aver")

#### we cluster and create a heatmap using heatmap.2 from gplots package
## we create a vector with colors that we will identify the treatment levels by
trt1 = c(rep("deepskyblue",5),rep("magenta",4))
## create a .pdf file
pdf(file="heatmap.pdf",height = 10,width=10)
## plot a heatmap and clustering using heatmap.2 using the transformed normalized data (t(norm.mt,[1:20])
## showing only the first 20 otus, show a dendrogram on the columns (now the samples) using the clustering
## we have already done above (Colv = as.dendrogram(row.clus)), color the sides of that dendrogram
## to reflect the tretment origin of the samples (ColSideColors=trt1), and put that dendrogram on the 
## column (dendrogram="col"), don't add a key and son't add a trace plot.
heatmap.2(t(norm.mt[,1:20]), Colv = as.dendrogram(row.clus), ColSideColors = trt1,dendrogram = "col",key = FALSE,trace="none")
## save
dev.off()

#### We cluster and plot a heatmap using pheatmap
## first we need our treatment levels that come from the meta.df data frame 
Trt = meta.df[,c(2,1)]
## create a .pdf file
pdf(file="pheatmap.pdf",width=10,height=10)
## use function pheatmap.3.ftn (from functions.R) to create a heatmap using the transformed normalzed matrix
## (t(norm.mt)) and treatment lvels found in Trt, assume number of clusters to be 3 (cutoff=3; and you can change this)
## cluster the samples using the above clustering (col.clust=row.clus), deside on cell height and width (cellh=3, cellw=10),
## and choose a font size (fonts = 7)
pheatmap.3.ftn(d.df = t(norm.mt),trt.df = Trt,cutoff = 3,col.clust = row.clus,cellh = 3,cellw = 10,fonts = 7)
## save
dev.off()

####### We can make the above plots based on taxonomic levels other than otus
#### to do so we need to create a phyloseq container based on the normalized data
#### Making sure that the normalized data is in a data frame format
data.df.nrm = data.frame(norm.mt)
#### choosing only the taxa observed in the normalized data to include in the phyloseq container
taxa.mt.nrm = as.matrix(taxa.df[colnames(data.df.nrm),])
#### the experimental data didn't change so we ad it as is
data.ps = phyloseq(otu_table(data.df.nrm,taxa_are_rows = FALSE),tax_table(taxa.mt.nrm),sample_data(meta.df))

#### Reduce the normalized data to the family level and extract the new data, taxa, and experimental design
data.ps.nrm.f = tax_glom(data.ps,"family")
data.df.nrm.f = data.frame(otu_table(data.ps.nrm.f))
taxa.df.nrm.f = data.frame(tax_table(data.ps.nrm.f)[,1:5])
colnames(data.df.nrm.f) = taxa.df.nrm.f$family
meta.df.nrm.f = data.frame(sample_data(data.ps.nrm.f))
norm.mt.f = as.matrix(data.df.nrm.f)

#### compute the ditance matrix and clustering (see above)
data.dist <- vegdist(norm.mt.f, method = "bray")
row.clus <- hclust(data.dist, "aver")

### draw the heatmap and clustering using pheatmap (see above) on the family level
Trt = meta.df.nrm.f[,c(2,1)]
pdf(file="pheatmap-family.pdf",width=5,height=5)
pheatmap.3.ftn(d.df = t(norm.mt.f),trt.df = Trt,cutoff = 3,col.clust = row.clus,cellh = 5,cellw = 10,fonts = 7)
dev.off()
