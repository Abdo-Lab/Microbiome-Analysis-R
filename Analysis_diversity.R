####### If you are running the Data processing step from here (source(file="Data_Processing.R) below)
#### then running the libraries.R (source(file="libraries.R")) is redundant but doesn't hurt
source(file="libraries.R")

######## Data Processing
#### If you are happy with your data processing step you can run it again here from the "Data_Processing.R" file without
#### having to open the file again
source(file = "Data_Processing.R")

#### We need row sums to calculate the rarefied richness below
d.rs = apply(data.df,1,sum)

######## Diversity analysis (old fashioned ecological approach [richness and diversity])
######## Done on the otu level!
#### We use estimate_richness() function in phyloseq to estimate diversity measures (here on non normalized data)
#### clickon d.rich to see what measures we have estimated
d.rich = estimate_richness(data.ps)
#### We use rarefy() function in vegan to normalize and calculate the expected richness based on the minimum
#### sequencing depth, and also compute the standard error se=TRUE (rarefy(data.df,min(d.rs),se=TRUE)) 
#### then we transform that calculated richness using t() (switching the rows and columns) 
rich.mt = t(rarefy(data.df,min(d.rs),se=TRUE))
#### We create a new data frame that incorporates columns 1, 6 and 8 from the richness estimates from phyloseq 
#### which are the observed richness, shannon diversity in InvSimpson. We add to that the first column (the expected richness)
#### from the vegan output and put it rich.df
rich.df = data.frame(cbind(d.rich[,c(1,6,8)],rich.mt[,1]))
#### We rename the columns to refelct what they are
names(rich.df) = c("Observed","Shannon","InvSimpson","Richness")

#### This step is in preparation to plot. We can't plot the above measures all in the same graph
#### unless we can have all the values in one column and partition by treatment level, and diversity measure.
#### We first create an empty container rich.plt = c()
rich.plt = c()
#### We extract the column names from the rich.df data frame and put them in rich.nm
rich.nm = colnames(rich.df)
#### We have for columns and hence 4 names so we iterate from 1 to 4 (length(rich.nm))
for(i in 1:length(rich.nm)){
  #### We creat a variable rich that includes the name of a diversity measure every iteration
  #### repeated as many as the sample number (same as length of column of rich.df [length(rich.df[,1])])
  rich = rep(rich.nm[i],length(rich.df[,1]))
  #### We combine the columns of the meta.df, the new rich and the diversity measure we using cbind() function
  #### We then add the combined matrix (data frame realy) to the end of the new data frame rich.plt using rbind
  rich.plt = rbind(rich.plt,cbind(meta.df,rich,rich.df[,i]))
}
#### We rename the new columns to reflect what they are (click on rich.plt after the next step to see what you have)
colnames(rich.plt) = c("sample","trt","rich","value")

#### We create a new .pdf file
pdf(file="diversity.pdf",width=5,height = 5)
#### qplot() is a function in ggplot2
#### We plot the diversity measures that are in rich.plt (data=rich.plt)
#### x axes is the treatment level, y axes is the diversity values and the type of plot is the geom and is a boxplot
#### Then we facet (split) the plots by diversity measure. 
qplot(data=rich.plt, x = trt, y = value, geom="boxplot")+facet_wrap(~rich,scale="free")
#### save
dev.off()

#### Testing richness using analysis of variance approaches
## First we combine the diversity data frame (rich.df) with the experimental design (meta.df) and we add the sequencin 
## depth to the mix (we add columns together using cbind())
depth = d.rs
rich.df1 = cbind(meta.df,rich.df,depth)

#### vegan richness
## We use lm (linear model) to fit an ANOVA model using treatement levels (A and B) that is we 
## fit the computed Ricness using vegan to fit a linear model (trt has to be factor)
## and this linear model is fit using the data in the new combined data frame rich.df1
rich.R.lm = lm(Richness~trt,data=rich.df1)
## We compute the ANOVA table based on the linear model fit
anova(rich.R.lm)
## We plot the standardized residuals to see how good the fit has been
## function residuals() finds the residuals, function sd() calculates the
## standard deviation for those residuals. And dividing the residuals 
## by their standard deviations results in stadardizing them
plot(residuals(rich.R.lm)/sd(residuals(rich.R.lm)))

#### Shannon
## Same as for Richness
rich.R.lm = lm(Shannon~trt,data=rich.df1)
anova(rich.R.lm)
plot(residuals(rich.R.lm)/sd(residuals(rich.R.lm)))

#### InvSimpson
## Same as for richness
rich.R.lm = lm(InvSimpson~trt,data=rich.df1)
anova(rich.R.lm)
plot(residuals(rich.R.lm)/sd(residuals(rich.R.lm)))

############# We can do the same after normalizing the data using the minimum sequencing depth
#### (phyloseq suggests doing the above without normalization)
#### I wouldn't use the subsampling approach
#### To normalize: 1) divide counts per sample by sequencing depth
p.df = data.df/d.rs
#### 2) create a new normalized count data frame by multiplying by minimum sequencing depth
d.norm.df = p.df * min(d.rs)
#### 3) create a new phyloseq container using the normalized data
d.norm.ps = phyloseq(otu_table(d.norm.df,taxa_are_rows = FALSE),as.matrix(taxa.df),sample_data(meta.df))
#### repeat all the above replacing data.ps with d.norm.ps and data.df with d.norm.df