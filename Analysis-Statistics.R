####### If you are running the Data processing step from here (source(file="Data_Processing.R) below)
#### then running the libraries.R (source(file="libraries.R")) is redundant but doesn't hurt
source(file="libraries.R")
#### We will use some of the functions I programmed for you (have a look in the functions.R file to see what they do)
source(file="functions.R")

######## Data Processing
#### If you are happy with your data processing step you can run it again here from the "Data_Processing.R" file without
#### having to open the file again
source(file = "Data_Processing.R")
#### adding the interaction to meta.df
meta.df$trt12 = paste(meta.df$trt1,meta.df$trt2,sep="")
#### switching columns to have sample id at end
meta.df = meta.df[,c(1,2,4,3)]

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

##### Redundancy Analysis
### Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible first order treatment effect using equation linking the 
### data matrix after standardizing to and the general term "." (the equation is hell.norm.mt ~ .), the data is represents the 
### two factors that we are using here (trt1 and trt2). we use the function rda() in vegan
mod1 <- rda(hell.norm.mt ~ trt1*trt2, data = meta.df[,1:2])
### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data = meta.df[,1:2])
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm")

### We look at the ANOVA table for the resulting best model
mod$anova
anova(mod)
### We evaluate significane using perMANOVA here again (calculate the p-values based on a permutation test)
anova(mod, by = "term", perm = 1000)

#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
ordiplot(mod,type="t")

##### distance-based Redundancy Analysis (dbRDA)
### We don't need to transform using hellinger or anything, other wise the steps are the same as above
### using the function dbrda() instead of rda() in vegan
mod1 <- dbrda(norm.mt ~ trt1*trt2, data = meta.df[,1:2],distance = "bray")
mod0 <- dbrda(hell.norm.mt ~ 1, data = meta.df[,1:2],distance = "bray")
mod <- step(mod0, scope = formula(mod1), test = "perm")

mod$anova
anova(mod)
anova(mod, by = "term", perm = 1000)


#######  metagenomeSeq analysis otu level
####  ZIG-Normal
#### Creating the model matrix (here it is a cell means model)
mod = model.matrix(~ -1 + trt12, B=1000, meta.df)
#### removing taxa that might be scarcely available
d.ms.fltr = filterData(d.ms,present = 10,depth = 1)
#### fitting the zero inflated gaussian model (ZIG) using the above specified model
d.fit = fitZig(d.ms.fltr, mod)
#### recovering the model parameters and significance tests
coef.df = MRcoefs(d.fit,number=94)
#### ordering the taxa based on their adjusted significance
coef.df = coef.df[order(coef.df$adjPvalues),]
#### Can output this table
write.csv(x,file="ZIN-otu.csv")
#### also can out put the names of those taxa
t.coef = taxa.df[row.names(coef.df[coef.df$adjPvalues<0.05,]),]
write.csv(t.coef,file="ZIN-otu-names.csv")
#### extracting model fit and design matrix to use to compare certain means
zigFit = d.fit$fit
finalMod = d.fit$fit$design
#### creating a contrast matrix to use to compare certain means
contrast.matrix1 = makeContrasts(trt12AC-trt12BC,trt12AD-trt12BD, levels = finalMod)
contrast.matrix2 = makeContrasts(((trt12AC+trt12AD)/2)-((trt12BC+trt12BD)/2), levels = finalMod)
c.fit1 = contrasts.fit(zigFit, contrast.matrix1)
c.fit1 = eBayes(c.fit1)
c.fit2 = contrasts.fit(zigFit, contrast.matrix2)
c.fit2 = eBayes(c.fit2)

#### Output of the contrast analysis (one of them)
coef.df.c2 = topTable(c.fit2, number=94)
write.csv(coef.df.c2,file="ZIN-otu-2.csv")
t.coef.df.c2 = taxa.df[row.names(coef.df.c2[coef.df.c2$adj.P.Val<0.1,]),]
write.csv(t.coef.df.c2,file="ZIN-otu-names-2.csv")

#### setup for plotting of significant contrasts
### level of significance
sig = 0.1
### chossing significant taxa
c.fit2.logFC = coef.df.c2$logFC[coef.df.c2$adj.P.Val<sig] 
c.fit2.genus = as.character(t.coef.df.c2$genus[coef.df.c2$adj.P.Val<sig])
x.col = c()
for(i in 1:length(c.fit2.logFC)){
  x.col = ifelse(c.fit2.logFC>0,"red","blue")
}
pdf(file="logFC.pdf",height=7.5,width=7.5)
par(mar=c(5,15,5,5))
barplot(c.fit2.logFC,names.arg = c.fit2.genus,col=x.col,horiz=TRUE,las=1,xlab = "log Fold Change",cex.names = 0.8)
dev.off()

