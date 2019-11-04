#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("Mfuzz")
rm(list=ls())
library('Mfuzz')
#read data
file<- read.table('fpkm.xls',header=T,row.names = 1)
#Set the missing value of the data to NA
file[file==0]=NA
#Mfuzz clustering requires an ExpressionSet type object, 
#so you need to build such an object first with expressions.
count<-data.matrix(file)
eset<-new("ExpressionSet",exprs=count)
#Remove genes that have more than one NA
eset<- filter.NA(eset, thres=0)
#Genes with small differences between samples were removed according to standard deviation
eset<-filter.std(eset,min.std=0)
#2standardized
eset<-standardise(eset)

#The number of clustering
c<-6
#Evaluate the best m value
m<-mestimate(eset)
#clustering
c1<-mfuzz(eset,c=c,m=m)

#visualization
mfuzz.plot(
  eset,
  c1,
  mfrow=c(2,3),
  new.window=FALSE)
