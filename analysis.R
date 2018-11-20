library(readxl)
library(metabolomics)
library(pca3d)
library(pcaMethods)
library(ggplot2)
library(rgl)
library(stringr)

# problem reading in data...variables are factors

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
forimput <- temp[,c(6689,1:6688)]
nomiss <- MissingValues(gooddata[,-c(1:2)],column.cutoff = 0.95,group.cutoff = 0.7)
# type summary(nomiss) shows 39 groups??

# summarize data
# need nonmissing data
sumG <- GroupSummary(nomiss)
##creating a vector log.mean to store the log mean of our QC samples
log.mean<-log(sumG$mean)

log.mean<-log(sumG$mean[4,])
#create a vector Rsd 
Rsd<-sumG$cv[4,]*100
#create a data frame vis.data with two columns  
vis.data<-data.frame(log.mean,Rsd)
#initialize ggplot object# set plot aesthetics
p<-ggplot(vis.data,aes(x=log.mean,y=Rsd))
# to plot points for scatterplot
p<-p+geom_point(alpha=.75,size=1)
print(p)

