setwd("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\")

library(readxl)
library(metabolomics)
library(pca3d)
library(pcaMethods)
library(ggplot2)
library(rgl)
library(stringr)
library(data.table)
library(dplyr)

# read and transpose data
mydata <- read_excel(path="H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\raw data for metaboanalyst no null.xlsx"
                     ,sheet=1 )
mydata <- as.data.frame(mydata)
mydata <- mydata[!is.na(mydata$Name),]
alldata <- as.data.frame(t(mydata[,-1]))
colnames(alldata) <- mydata$Name

# create variable for ID
for (i in 1:nrow(alldata)) {
  alldata$longid[i] <- row.names(alldata)[i]
}
alldata$longid <- gsub("\\s", "", alldata$longid) 
alldata$id <- gsub("OGTT0", "", alldata$longid) 
alldata$id <- gsub("BCTP0", "", alldata$id) 
#row.names(alldata) <- alldata$id
# not working to set rownames equal to ID b/c doesn't allow duplicate rownames
# how else do you specify paired nature of data?
View(alldata[c("longid","id")])

# convert to numeric
num <- alldata[ ,!(colnames(alldata) %in% c("Batch","id","longid"))]
fwrite(num,"C:\\Temp\\temp.csv")
newnum <- fread("C:\\Temp\\temp.csv",colClasses = "numeric")
alldata <- alldata[ ,(colnames(alldata) %in% c("Batch","id"))]
alldata <- cbind(alldata,newnum)

# find variables that have too few nonmissing observations
check <- as.data.frame(apply(alldata,2,function(x) length(which(!is.na(x)))))
for (i in 1:nrow(check)) {
  check$id[i] <- row.names(check)[i]
}
bad <- as.data.frame(check[check[,1]<3,])
gooddata <- alldata[ ,!(colnames(alldata) %in% bad$id)]
gooddata$Group[gooddata$Batch=="diet"] <- 1
gooddata$Group[gooddata$Batch=="No-diet"] <- 0
for (i in 1:nrow(gooddata)) {
  row.names(gooddata[i,]) <- gooddata$id[i]
}
# also get rid of 3 participants who only have one timepoint
gooddata <- gooddata[gooddata$id !='Control6164-23' & gooddata$id !="PCOS6164-28" & 
                       gooddata$id !="PCOS6164-42",]
# sort so that participants are in the same order within condition
gooddata <- arrange(gooddata,Group,id)

# create a dataset using only the good data but in the format for metabolomics package
temp <- gooddata[,-c(1:2)]
cn <- ncol(temp)
gooddata.format <- temp[,c(cn,1:(cn-1))]

# impute missing data
nomiss <- MissingValues(gooddata.format,column.cutoff = 0.95,group.cutoff = 0.7,saveoutput = TRUE,
                        outputname = "C:\\Temp\\newoutput",complete.matrix = TRUE)
nomissdf <- fread("C:\\Temp\\newoutput.csv",header=TRUE)
for (i in 1:nrow(nomissdf)) {
   row.names(nomissdf)[i] <- nomissdf[i,1] 
}

# log transform gooddata
gooddata.log <- LogTransform(gooddata.format)$output

#Separating by diet/no diet
gooddata.group<-factor(gooddata.log[,1],levels=unique(gooddata.log[,1]))
dietmat<-gooddata.log[which(gooddata.log[,1]==1),-1]
nodietmat<-gooddata.log[which(gooddata.log[,1]==0),-1]

#Linear model fit with ordinary statistics
ordFit<-LinearModelFit(datamat=data.matrix(dietmat-nodietmat),
                       ruv2=FALSE,
                       factormat=matrix(1,nrow=nrow(dietmat)))
TwoGroupPlots(gooddata.log[,-1],
              tstats = ordFit$t[,1],
              foldchanges = ordFit$coef[,1],
              pvalues = ordFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)


#Linear model fit with moderated statistics
modFit<-LinearModelFit(datamat=data.matrix(dietmat-nodietmat),
                       ruv2=FALSE,
                       moderated=TRUE,
                       factormat=matrix(1,nrow=nrow(dietmat)))
TwoGroupPlots(gooddata.log[,-1],
              tstats = modFit$t[,1],
              foldchanges = modFit$coef[,1],
              pvalues = modFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)

# HeatMap - not working
HeatMap(nomissdf,colramp=redgreen(75),scale="column",dendrogram = "both",distmethod = "euclidean",
        aggmethod = "complete",key=TRUE,margins=c(1,2),ColSideColors = "Group")

# prep the data # data pre-treatement autoscale
md<-prep(log(gooddata.format[,-c(1)]),scale="uv",center=T)
# prcomp function for principal components - doesn't work b/c of missing data
# probject<-prcomp(~.,data=md,na.action=na.pass,scale=TRUE)
# need to use NIPALS PCA due to missing data
a <- checkData(as.matrix(gooddata.format))
probject <- pca(md,method="nipals",nPcs = 3)
plotPcs(probject,pcs=1:3,type="scores",col=as.factor(gooddata.format$Group))

# need to figure out how to do PLS-DA