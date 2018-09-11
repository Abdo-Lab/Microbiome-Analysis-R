##### Load all needed packages that will be used in all in analyses that you 
##### can find in the file libraries.R (the function "source" reads and 
##### executes all commands in another file). Libraries have to be loaded 
##### every time you run the following programs (they don't need to be installed; 
##### they are installed only once after installing or updating R)
source("libraries.R")

############# Uploading mothur's results
#### read.table: reads a tab delimited file. Here is is reading the samples.shared file resulting from 
#### mothur that is located in the directory "results". The file includes rows representing the samples
#### and columns that represent the otus. We use the otu names to be the column names by choosing "header = TRUE"
#### click on data.df (right upper pane) to see how it looks like
data.df = read.table(file="results/samples.shared",header=TRUE)

#### loading taxonomy using phyloseq in different columns 
#### We use the command "import_mothur" with option "mothur_costaxonomy_file = " to read and store 
#### the mothur taxonomic assignment available in file "samples.0.03.cons.taxonomy" also saved in "results" 
#### this results in a phyloseq container (taxa.phs) and splits the taxonomic assignment into the different taxonomic levels
#### and add them each in a separate column. It will allow us to do analyses on the taxonomic level we choose using some 
#### manipulations that are first described in Analysis_exploratory.R
taxa.phs = import_mothur(mothur_constaxonomy_file = "results/samples.0.03.cons.taxonomy")
#### The following uses colnames to rename the columns of the taxonomic table created above to refelct the taxonomic levels
colnames(tax_table(taxa.phs)) = c("kingdom", "phylum", "class", "order", "family",  "genus")
#### The following changes the format of the taxonomic table from a phyloseq container to a data.frame
#### it first uses tax_table() a function in phyloseq to extract the table and then uses data.frame (or as.data.frame)
#### to change the format to a data frame
#### click on the taxa.df (right upper pane) to see how it looks like
#### the otu names form the row names of this data frame 
taxa.df = data.frame(tax_table(taxa.phs))

############## manipulatin the "results" files
#### We save the sample names (column # 2) from the current data.df data frame
r.nm = as.character(data.df[,2])
#### We use data.frame() function to replace data.df by a new data frame that removes the 
#### first 3 columns (data.df[.-c(1:3)]) and that uses the sample names r.nm as row names (row.names=r.nm)
data.df = data.frame(data.df[,-c(1:3)],row.names=r.nm)
#### We take the column names (otu00001, otu00002, ...) and save them into otu.nm
otu.nm = names(data.df)

#### The following makes sure that the order of the rows of the taxa.df (the taxonomy table) 
#### exactly matches the order of the columns of the otu table using row.names(taxa.df) to extract the row names
#### of taxa.df and comparing them using %in% to otu.nm that includes the column names from the data.df
taxa.df = taxa.df[row.names(taxa.df)%in%otu.nm,]

#### the following is only useful if you have a Mock community and that you named the sample that includes it "Mock"
#### If you have a mock community with that you named something else (like ZymoSTD for example) you should change "Mock"
#### to that name ("ZymoSTD" for example)
## Here we choose the row from data.df that belongs to the "Mock" community. We first look at all data.df row names 
## using row.names(data.df) and we use == to find the exact match to the word "Mock" (row.names(data.df)=="Mock").
## This comparison is put in the first index of the data.df[(row.names(data.df)=="Mock",]
## At the same time we choose to only look at the columns (otus) that have counts larger than zero part of the
## "Mock" community so we choose the "Mock" row from data.df using "data.df[row.names(data.df)=="Mock",]" and compare 
## that whole row to 0 to find only the values greater than 0 (data.df[row.names(data.df)=="Mock",]>0) and add this to 
## the second index data.df[row.names(data.df)=="Mock",data.df[row.names(data.df)=="Mock",]>0]
data.df[row.names(data.df)=="Mock",data.df[row.names(data.df)=="Mock",]>0]

#### Using what we see from the above command we choose a cutoff (a number that we use to subtract from all counts)
#### Rational for this number is as follows: We know that we have about 21 bacterial species (8 using the Zymo standard)
#### that form the Mock community so it doesn't make sense to see more (or less) than that number of otus in the 
#### Mock samples. That meeans that if the total number is larger than 20-21 (or larger than 8 using Zymo standard) then 
#### we have artifacts due to sequencing error and/or contamination and we choose a cutoff that will minimize that effect.
#### Here we chose 1 as our cutoff which reduces the total number of taxa we observe in the Mock to about 19.
## 1) subtract cutoff zeros will become -1, 1's will become 0's 2's will become 1's and so on and put the outcome
##    in a new data frame (d.df)
d.df = data.df - 2
## 2) change all negative values in the table to 0's (find those negative values by searching the whole data frame
##    for them using the condition d.df<0)
d.df[d.df<0] = 0
## 3) find the column sums; we want to find those otus that disappeared after this new error correction
##    (like we said in the primer to R lecture function "apply" will apply the "sum" to the columns "# 2" 
##    of the data frame "d.df"; apply(d.df,2,sum))  
d.cs = apply(d.df,2,sum)
## 4) keep all columns (otus) that have sums greater than 0
d.df = d.df[,d.cs > 0]
## 5) keep only the otus that match the otu table (d.df) by choosing these rows with the same names as the column names 
##    of the otu (d.df) table and hence removing all otus that have 0 abundance accross all samples 
t.df = taxa.df[colnames(d.df),]

############ Metadata
#### to do any data analysis we need an experimental design and that comes from the "samples.files" file that we adjust
#### either using R as below or using excel or any text editor
#### first we read this file into "meta.df" using read.tables (this file doesn't have a header or column names; it only has a sample names 
#### and file names). We read it as is (using as.is=TRUE) just to make sure it is loaded as text an not converted to factors 
#### or numbers
meta.df = read.table(file="samples.files",header=FALSE,as.is = TRUE)
#### The first column is the sample names. We use the function "data.frame" to make sure that meta.df is a data frame with row names 
#### that equal to the sample names (these should exactly match the row names of the transfomred d.df)
#### remember columns of d.df (and data.df) is linked to rows t.df (and taxa.df) by otu names and rows of d.df (and data.df) and 
#### rows of meta.df are linked by the sample names
meta.df = data.frame(meta.df[,1],row.names=meta.df[,1])
#### life is made easier if both meta.df and d.df have the same sample-name order. First column still has sample names 
#### (row names are also sample names)
d.df = d.df[row.names(meta.df),]

#### Now we are ready to remove the Mock community data and we will not use Mock or any csequencing control after this point
#### In this data (click on d.df to see) Mock couns are in row 10 of the d.df data frame and we remove that row by choosing
#### only the lines that do not include line 10 (-10 in the first index). We take these lines away from both the d.df and meta.df
#### we match by row
d.df = d.df[row.names(d.df)!="Mock",]
meta.df = meta.df[row.names(d.df),,drop=FALSE]
#### Now create the experimental design (can do this in excel)
#### We have two treatments A (for eaxmple "Sick") and B (for example "Control")
#### We have 9 samples we claim that 5 belong to level A (sick) and 4 belong to level B (Control)
#### we create a new variable "trt" that represents treatment "trt"
trt = c(rep("A",5),rep("B",4))
#### we make it be a factor using as.factor (or factor)
trt = as.factor(trt)
#### and add it in a new column (the last column) that is called trt (click on meta.df after the next step to see what happened)
meta.df$trt = trt
#### We rename the columns to "trt" and "id" (sample id)
colnames(meta.df) = c("trt","id")

#### Now data is clean and well organized and we replace data.df with d.df and taxa.df with t.df (we didn't do this earlier 
#### to avoid loading the data incase we made a mistake)
data.df = d.df
taxa.df = t.df

########### We use the data.df (the otu table), taxa.df (the taxa table) and the experimental design meta.df (sample data)
########### to create a new phyloseq container
### 1) we use otu_table to incorportate data.df telling phyloseq that the otus are the columns and not the rows,
### 2) we temporary change taxa.df from a data.frame to a matrix and incorporate it as the taxonomic table (tax_table(as.matrix(taxa.df)))
### (we can also use tax_table(taxa.phs) to incorporate thse taxa)
### 3) we use sample_data() to incorporate the experimental design
data.ps = phyloseq(otu_table(data.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.df)),sample_data(meta.df))

#### All what we need to keep at this point are the data.df, taxa.df and meta.df data frames along with the 
#### phyloseq container data.ps, so we remove the rest
rm(d.cs,r.nm,otu.nm,d.df,t.df,taxa.phs,trt)
