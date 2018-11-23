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

# summarize data
# this section of code has to do with QC samples, and don't have those in my file
# so I modified this to use the actual samples
sumG <- GroupSummary(nomissdf[,-1])
##creating a vector log.mean to store the log mean of our samples for each group
log.meangp0 <- log(sumG$means[1,])
log.meangp1 <- log(sumG$means[2,])
#create a vector Rsd 
Rsdgp0 <-sumG$cv[1,]*100
Rsdgp1 <-sumG$cv[2,]*100
#create a data frame vis.data with two columns  
vis.datagp0 <- data.frame(log.meangp0,Rsdgp0)
vis.datagp1 <- data.frame(log.meangp1,Rsdgp1)
#initialize ggplot object# set plot aesthetics
p0 <-ggplot(vis.datagp0,aes(x=log.meangp0,y=Rsdgp0))
# to plot points for scatterplot
p0<-p0+geom_point(alpha=.75,size=1)
print(p0)
#initialize ggplot object# set plot aesthetics
p1 <-ggplot(vis.datagp1,aes(x=log.meangp1,y=Rsdgp1))
# to plot points for scatterplot
p1<-p1+geom_point(alpha=.75,size=1)
print(p1)

# example to show Normalize function
mydata.log <- LogTransform(nomissdf[,-1])$output
for (i in 1:nrow(mydata.log)) {
  row.names(mydata.log)[i] <- nomissdf[i,1] 
}
norm_is <- Normalise(mydata.log, method = "median",saveoutput = TRUE,outputname = "C:\\Temp\\norm")
norm <- fread("C:\\Temp\\norm.csv",header=TRUE)
for (i in 1:nrow(nomissdf)) {
  row.names(norm)[i] <- nomissdf[i,1] 
}
# now summarize data after normalization
# summarize data
# this section of code has to do with QC samples, and don't have those in my file
# so I modified this to use the actual samples
sumG_norm <- GroupSummary(norm[,-1])
##creating a vector log.mean to store the log mean of our samples for each group
meangp0_norm <- (sumG_norm$means[1,])
meangp1_norm <- (sumG_norm$means[2,])
#create a vector Rsd 
Rsdgp0_norm <-sumG_norm$cv[1,]*100
Rsdgp1_norm <-sumG_norm$cv[2,]*100
#create a data frame vis.data with two columns  
vis.datagp0_norm <- data.frame(meangp0_norm,Rsdgp0_norm)
vis.datagp1_norm <- data.frame(meangp1_norm,Rsdgp1_norm)
#initialize ggplot object# set plot aesthetics
p0_norm <-ggplot(vis.datagp0_norm,aes(x=meangp0_norm,y=Rsdgp0_norm))
# to plot points for scatterplot
p0_norm <- p0_norm+geom_point(alpha=.75,size=1)
print(p0_norm)
#initialize ggplot object# set plot aesthetics
p1_norm <-ggplot(vis.datagp1_norm,aes(x=meangp1_norm,y=Rsdgp1_norm))
# to plot points for scatterplot
p1_norm <- p1_norm+geom_point(alpha=.75,size=1)
print(p1_norm)
# normalized data plots look very strange
# is it possible to normalize using internal or external standards as done in demo?

# paired comparison of non-normalized data
comp <- TwoGroup(mydata.log, alternative = "two.sided",paired=TRUE)
# think the IDs are wrong - need to be the same across conditions