# download the BioC installation routines
source("http://bioconductor.org/biocLite.R")
# install the core packages
biocLite()
# install the GEO libraries
biocLite("GEOquery")
biocLite("makecdfenv")
biocLite("simpleaffy")
biocLite("hgu133a.db", type = "source")

install.packages("tibble")
library(affy)
library(makecdfenv)
library(limma)

library(GEOquery)
library(simpleaffy)

library(hgu133a.db)
library(hgu133acdf)
setwd("data")
data = ReadAffy(cdfname = "hgu133acdf")

data_norm = rma(data)

rma1 = exprs(data_norm)

rma1[1:5, 1:5]


# Write RMA-normalized, mapped data to file
write.table(rma, file = "~/R Project/rma.txt", quote = FALSE, sep = "\t")

human4302mmug<-make.cdf.env("HGU133A_Hs_ENSG.cdf")


setwd("data")
#category<-read.table("classify.txt",header=TRUE,fileEncoding="utf-16",sep= "\t")

category<-read.table("classify.txt",header=TRUE)
rawdata<-ReadAffy(phenoData = category,cdfname = "human4302mmug")

eset<- rma(rawdata)
data_matrix1<-exprs(eset)

boxplot(exprs(eset))
write.exprs(eset,file="~/R Project/classifiexpr.txt")


getwd()
#diffrential expression analysis
eset$TYPE
#subsetting our data

ourdata<-eset[,eset$TYPE %in% c("Normal","Osteoarthritis","Rheumatoid_arthritis")]

ourdata
ourdata$TYPE

ourdata$TYPE<-factor(ourdata$TYPE)

#desingn
design2<-model.matrix(~ ourdata$TYPE-1)
head(design2)
colnames(design2)<-c("Normal","Osteoarthritis","Rheumatoid_arthritis")


fit3<-lmFit(ourdata,design2)
contrast<-makeContrasts("Normal-Osteoarthritis-Rheumatoid_arthritis",levels = design2)

fit4<-contrasts.fit(fit3,contrast)
fit5<-eBayes(fit4)


#viewing the results for top 50 differentially expressed genes you can change
topTable(fit5,number=50)