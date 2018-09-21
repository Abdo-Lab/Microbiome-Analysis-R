##### Load all needed packages that will be used in all in analyses that you 
source("libraries.R")

############# Uploading mothur's results
data.df = read.table(file="results/samples.shared",header=TRUE)

#### loading taxonomy using phyloseq in different columns 
taxa.phs = import_mothur(mothur_constaxonomy_file = "results/samples.0.03.cons.taxonomy")
colnames(tax_table(taxa.phs)) = c("kingdom", "phylum", "class", "order", "family",  "genus")
taxa.df = data.frame(tax_table(taxa.phs))

############## manipulatin the "results" files
#### We save the sample names (column # 2) from the current data.df data frame
r.nm = as.character(data.df[,2])
data.df = data.frame(data.df[,-c(1:3)],row.names=r.nm)
otu.nm = names(data.df)

#### The following makes sure that the order of the rows of the taxa.df (the taxonomy table) 
taxa.df = taxa.df[row.names(taxa.df)%in%otu.nm,]

#### the following is only useful if you have a Mock community and that you named the sample that includes it "Mock"
#### If you have a mock community with that you named something else (like ZymoSTD for example) you should change "Mock"
#### to that name ("ZymoSTD" for example)

a = data.df[row.names(data.df)=="Mock",data.df[row.names(data.df)=="Mock",]>0]
length(a)



#### Using what we see from the above command we choose a cutoff (a number that we use to subtract from all counts)
d.df = data.df - 1
d.df[d.df<0] = 0
d.cs = apply(d.df,2,sum)
d.df = d.df[,d.cs > 0]
t.df = taxa.df[colnames(d.df),]

############ Metadata
meta.df = read.table(file="sample.files",header=FALSE,as.is = TRUE)
meta.df = data.frame(meta.df[,1],row.names=meta.df[,1])
d.df = d.df[row.names(meta.df),]

#### Now we are ready to remove the Mock community data and we will not use Mock or any csequencing control after this point
d.df = d.df[row.names(d.df)!="Mock",]
meta.df = meta.df[row.names(d.df),,drop=FALSE]
#### Now create the experimental design (can do this in excel)
#### We have two treatments A (for eaxmple "Sick") and B (for example "Control")
trt = c(rep("A",5),rep("B",4))
#### we make it be a factor using as.factor (or factor)
trt = as.factor(trt)
#### and add it in a new column (the last column) that is called trt (click on meta.df after the next step to see what happened)
meta.df$trt = trt
#### We rename the columns to "trt" and "id" (sample id)
colnames(meta.df) = c("id","trt")

#### Now data is clean and well organized and we replace data.df with d.df and taxa.df with t.df (we didn't do this earlier 
#### to avoid loading the data incase we made a mistake)
data.df = d.df
taxa.df = t.df

########### We use the data.df (the otu table), taxa.df (the taxa table) and the experimental design meta.df (sample data)
########### to create a new phyloseq container
data.ps = phyloseq(otu_table(data.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.df)),sample_data(meta.df))

#### All what we need to keep at this point are the data.df, taxa.df and meta.df data frames along with the 
#### phyloseq container data.ps, so we remove the rest
rm(d.cs,r.nm,otu.nm,d.df,t.df,taxa.phs,trt,a)
#save.image(file="data.RData")
