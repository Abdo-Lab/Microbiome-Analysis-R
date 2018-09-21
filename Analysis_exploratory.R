####### If you are running the Data processing step from here (source(file="Data_Processing.R) below)
#### then running the libraries.R (source(file="libraries.R")) is redundant but doesn't hurt
source(file="libraries.R")
#### The following loads some functions we use to plot and look at our data (description is here and in the "functions.R" file)
source(file="functions.R")

######## Data Processing
#### If you are happy with your data processing step you can run it again here from the "Data_Processing.R" file without
#### having to open the file again
source(file = "Data_Processing.R")

######## Sequence depth
#### To assess sequencing depth per sample we have to sum over the rows of data.df (remember each row represents a sample
#### with row names equavilent to sample names)
#### we use apply to sum obtain row sums (# 1) from the data.df
rs.vt = apply(data.df,1,sum) # row sum (total reads per sample)

#### We creat a .pdf file and name it "depth.pdf" with dimentions 4 x 4 in
pdf(file="depth.pdf",width=4,height=4)
#### We plot a histogram representing the row sums "rs.vt" with main title = "Distribution of Sequence Depth and
#### x label = sequencing depth and y label = Number of Sample per that depth
#### this will be printed to the "depth.pdf" file
hist(rs.vt,main = "Distribution of Sequencing Depth",xlab = "sequencing depth",ylab = "Number of Samples")
#### We save the file
dev.off()
#### we also calculate the minimum and maximum sequencing depth to see how spread it is and assess need to normalize
min(rs.vt)
max(rs.vt)

######## Rarifaction curves (we discussed those in class)
##### We creat a .pdf file called rarefaction_curves.pdf with dimentions 10x10in (that you can change depending 
##### on the size you want and the resolution you have) to save the rarefaction curves
pdf(file="rarefaction_curves.pdf", width = 10, height = 10)
##### rarecurve() is a function in vegan it takes the otu table x=data.df, and to have different colors per samples
##### we use col and set these colors to the numbers from 1 to the total number of samples "col=1:length(rs.vt)"
##### label of the x axes is "Number of reads" and that for the y axes is "OTU"
rarecurve(x=data.df,col=1:length(rs.vt),xlab = "Number of reads",ylab="OTU")
#### We save the file
dev.off()

########### A look at proportions before normalizing per otu
####### The data frames created above were per otu and you might want to look at the distribution of these otus per experiment 
####### or per treatment level
## All samples bar plot (you already know how to create a .pdf file [see above])
pdf(file="barplots-allsamples.otu.pdf",width=10,height=10)
## We use my function bar.taxa.sample.ftn() (full description in "functions.R") to plot the barplot using the otu table data.df 
## We avoid plotting rare taxa by setting up a cutoff of 0.01 where anything that has a proportion of less than 1% per sample is
## left out of the plot (makes things go faster)
bar.taxa.sample.ftn(d.df=data.df,cutoff = 0.01)
## save
dev.off()

## Per treatment barplot for all otus
## create a .pdf file
pdf(file="barplots-per-treatment-otu.pdf",width=10,height=10)
## use bar.sample.trt.ftn() (full description in "functions.R") to provide bar plot per treatment level (A and B here)
## using data.df, and the experimental design (first column here must be treatment level and last column must be sample id)
## we also use a per sample cutoff of 1%. We call the figure "Barplot", the x label "Sample" and the y label "OTU"
bar.sample.trt.ftn(data=data.df,Trt=meta.df[,c(2,1)],cutoff = 0.01, title = "Barplot",xlab = "Sample",llab = "OTU")
## save file
dev.off()

####### It is cleaner to look at the barplots using a level higher than the otu level
####### we can choose to aollapse the otu tables (finding the counts per genus, family, order ... etc) using the 
####### tax_glom() function in phyloseq and this is why we created the phyloseq container data.ps
#### Here we apply the function tax_glom to data.ps to create a new phyloseq container that collapses the data to the 
#### family level (the otu table and the taxa table are now on the family level) we call the new container data.ps.f
data.ps.f = tax_glom(data.ps,"family")
#### We extract the new collapsed otu table and put it into a new data frame data.df.f using function otu_table in phyloseq 
data.df.f = data.frame(otu_table(data.ps.f))
#### We extract the new taxonomy table and save it in a new data frame taxa.df.f using function tax_table in phyloseq
taxa.df.f = data.frame(tax_table(data.ps.f)[,1:5])
#### We rename the columns of the data.df.f from otu00001, otu00004, ... to the family names
colnames(data.df.f) = taxa.df.f$family
#### We also extract the experimental design (which will not change) into data.df.f using function sample_data in phyloseq
meta.df.f = data.frame(sample_data(data.ps.f))

## Create a barplot for all samples on the family level (change llab to "Family" [default is "OTU"])
## and use per sample cutoff of waht to plot of 1% proportion (families with proportion above 1% will be included in the plot per
## sample).
pdf(file="barplots-allsamples.family.pdf",width=10,height=10)
bar.taxa.sample.ftn(d.df=data.df.f,cutoff = 0.01,llab = "Family")
dev.off()

## Create a barplot per treatment (A and B here) at the family level (see above)
pdf(file="barplots-per-treatment-family.pdf",width=5,height=5)
bar.sample.trt.ftn(data=data.df.f,Trt=meta.df.f[,c(2,1)],cutoff = 0.01, title = "Barplot",xlab = "Sample",llab = "Family")
dev.off()
