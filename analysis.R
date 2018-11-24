setwd("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\")

library(readxl)
library(metabolomics)
library(pca3d)
library(pcaMethods)
library(ggplot2)
library(rgl)
library(stringr)
library(data.table)

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
row.names(alldata) <- alldata$id
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
bad <- as.data.frame(check[check[,1]<9,])
gooddata <- alldata[ ,!(colnames(alldata) %in% bad$id)]
gooddata$Group[gooddata$Batch=="diet"] <- 1
gooddata$Group[gooddata$Batch=="No-diet"] <- 0
for (i in 1:nrow(gooddata)) {
  row.names(gooddata[i,]) <- gooddata$id[i]
}

# prep the data # data pre-treatement autoscale
md<-prep(log(gooddata[,-c(1:2)]),scale="uv",center=T)
# prcomp function for principal components - doesn't work b/c of missing data
# probject<-prcomp(~.,data=md,na.action=na.pass,scale=TRUE)
# need to use NIPALS PCA due to missing data
a <- checkData(as.matrix(gooddata))
probject <- pca(md,method="nipals",nPcs = 3)
plotPcs(probject,pcs=1:3,type="scores",col=gooddata$Group)

# impute missing data
temp <- gooddata[,-c(1:2)]
forimput <- temp[,c(6578,1:6577)]
nomiss <- MissingValues(forimput,column.cutoff = 0.95,group.cutoff = 0.7,saveoutput = TRUE,
                        outputname = "C:\\Temp\\newoutput",complete.matrix = TRUE)
nomissdf <- fread("C:\\Temp\\newoutput.csv",header=TRUE)
for (i in 1:nrow(nomissdf)) {
   row.names(nomissdf)[i] <- nomissdf[i,1] 
}

# paired comparison of non-normalized data
comp <- TwoGroup(mydata.log, alternative = "two.sided",paired=TRUE)
# think the IDs are wrong - need to be the same across conditions