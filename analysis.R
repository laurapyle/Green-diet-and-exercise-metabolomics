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
library(mixOmics)

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
alldata$id <- gsub("PCOS6164-", "", alldata$id)
alldata$id <- gsub("PCOSHS6164-", "", alldata$id)
alldata$id <- gsub("Control6164-", "", alldata$id)
alldata$uniqueid <- paste(substr(alldata$longid,1,1),alldata$Batch,alldata$id,sep="")
row.names(alldata) <- alldata$uniqueid
# not working to set rownames equal to ID b/c doesn't allow duplicate rownames
# how else do you specify paired nature of data?
View(alldata[c("Batch","longid","id","uniqueid")])

# convert to numeric
num <- alldata[ ,!(colnames(alldata) %in% c("Batch","id","longid","uniqueid"))]
fwrite(num,"C:\\Temp\\temp.csv")
newnum <- fread("C:\\Temp\\temp.csv",colClasses = "numeric")
alldata <- alldata[ ,(colnames(alldata) %in% c("Batch","id","uniqueid"))]
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
# also get rid of 3 participants who only have one timepoint
gooddata <- gooddata[gooddata$id !='23' & gooddata$id !="28" & 
                       gooddata$id !="42",]
# sort so that participants are in the same order within condition
gooddata <- gooddata[with(gooddata, order(Group, id)),  ]

# create a dataset using only the good data but in the format for metabolomics package
temp <- gooddata[,-c(1:3)]
cn <- ncol(temp)
gooddata.format <- temp[,c(cn,1:(cn-1))]

# impute missing data
nomiss <- MissingValues(gooddata.format,column.cutoff = 0.95,group.cutoff = 0.7,saveoutput = TRUE,
                        outputname = "C:\\Temp\\newoutput",complete.matrix = TRUE)
nomissdf <- fread("C:\\Temp\\newoutput.csv",header=TRUE)
nomissdf <- nomissdf[,-1]
for (i in 1:nrow(nomissdf)) {
  row.names(nomissdf)[i] <- row.names(gooddata.format)[i] 
}
nomissdf.log <- LogTransform(nomissdf)$output

# log transform gooddata
gooddata.log <- LogTransform(gooddata.format)$output

#Separating by diet/no diet
gooddata.group<-factor(gooddata.log[,1],levels=unique(gooddata.log[,1]))
dietmat<-gooddata.log[which(gooddata.log[,1]==1),-1]
nodietmat<-gooddata.log[which(gooddata.log[,1]==0),-1]

#Linear model fit with ordinary statistics
ordFit<-LinearModelFit(datamat=data.matrix(dietmat-nodietmat),
                       ruv2=FALSE,
                       factormat=matrix(1,nrow=nrow(dietmat)),outputname = "C:\\Temp\\ordFit",
                       saveoutput = TRUE)
TwoGroupPlots(gooddata.log[,-1],
              tstats = ordFit$t[,1],
              foldchanges = ordFit$coef[,1],
              pvalues = ordFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)
# I don't understand why some compounds have NA's for p-values

#Linear model fit with moderated statistics
modFit<-LinearModelFit(datamat=data.matrix(dietmat-nodietmat),
                       ruv2=FALSE,
                       moderated=TRUE,
                       factormat=matrix(1,nrow=nrow(dietmat)),outputname = "C:\\Temp\\modFit",
                       saveoutput = TRUE)
TwoGroupPlots(gooddata.log[,-1],
              tstats = modFit$t[,1],
              foldchanges = modFit$coef[,1],
              pvalues = modFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)
# trying TwoGroupPlot with anything significantly different, regardless of fold change
TwoGroupPlots(gooddata.log[,-1],
              tstats = modFit$t[,1],
              foldchanges = modFit$coef[,1],
              pvalues = modFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(1),
              pcutoff = 0.05)

# Volcano plot
VolcanoPlot(folds=modFit$coef, pvals=modFit$p.value,plimit=0.05)
# VolcanoPlot doesn't like NA's for p-values

# Box plots of specific metabolites that have NA for p-value
# Taurine is only present in one non-diet sample
# if not present in many samples in either group, not interesting - should we filter these out?
# if present in many samples in one group and few in the other, is this interesting?
MetBoxPlots(gooddata.log,"TAURINE (M-H)-")
MetBoxPlots(gooddata.log,"3-OXALOMALIC ACID (M-H)-[-H2O]")
MetBoxPlots(gooddata.log,"N-AMIDINO-ASPARTIC ACID (M+Cl)-")
MetBoxPlots(gooddata.log,"GLUCONIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"2,3-BISPHOSPHO-D-GLYCERIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"INOSINE 5'-DIPHOSPHATE (M-H)-")
MetBoxPlots(gooddata.log,"XANTHOSINE 5'-PHOSPHATE (2M-H)+")
MetBoxPlots(gooddata.log,"MALEIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"ORNITHINE (M-H)-")
MetBoxPlots(gooddata.log,"P-ACETAMIDOPHENYL BETA-D-GLUCURONIDE (M-H)-")
MetBoxPlots(gooddata.log,"1,7-DIMETHYL URIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"PHENYLPYRUVIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"PIMELIC ACID (M-H)-")
tapply(gooddata$"TAURINE (M-H)-",gooddata$Group, summary)
tapply(gooddata$"3-OXALOMALIC ACID (M-H)-[-H2O]",gooddata$Group, summary)
tapply(gooddata$"N-AMIDINO-ASPARTIC ACID (M+Cl)-",gooddata$Group, summary)
tapply(gooddata$"GLUCONIC ACID (M-H)-",gooddata$Group, summary)
tapply(gooddata$"2,3-BISPHOSPHO-D-GLYCERIC ACID (M-H)-",gooddata$Group, summary)
tapply(gooddata$"INOSINE 5'-DIPHOSPHATE (M-H)-",gooddata$Group, summary)
tapply(gooddata$"XANTHOSINE 5'-PHOSPHATE (2M-H)+",gooddata$Group, summary)
tapply(gooddata$"MALEIC ACID (M-H)-",gooddata$Group, summary)
tapply(gooddata$"ORNITHINE (M-H)-",gooddata$Group, summary)
tapply(gooddata$"P-ACETAMIDOPHENYL BETA-D-GLUCURONIDE (M-H)-",gooddata$Group, summary)
tapply(gooddata$"1,7-DIMETHYL URIC ACID (M-H)-",gooddata$Group, summary)
tapply(gooddata$"PHENYLPYRUVIC ACID (M-H)-",gooddata$Group, summary)
tapply(gooddata$"PIMELIC ACID (M-H)-",gooddata$Group, summary)

# Dendrogram
Dendrogram(gooddata.log)
# is there an interaction between diet and PCOS?  PCO group seem like they are more similar regardless of diet

# HeatMap 
HeatMap(nomissdf,colramp=redgreen(75),margins = c(5,10),key=FALSE)

# PCA plot
PcaPlots(nomissdf.log,scale=TRUE, center=TRUE)

# RLA plot
RlaPlots(gooddata.log,type="ag")
RlaPlots(gooddata.log,type="wg")

# prep the data # data pre-treatement autoscale
md<-prep(log(gooddata.format[,-c(1)]),scale="uv",center=T)
# prcomp function for principal components - doesn't work b/c of missing data
# probject<-prcomp(~.,data=md,na.action=na.pass,scale=TRUE)
# need to use NIPALS PCA due to missing data
a <- checkData(as.matrix(gooddata.format))
probject <- pca(md,method="nipals",nPcs = 3)
plotPcs(probject,pcs=1:3,type="scores",col=as.factor(gooddata.format$Group))

# PLS-DA
splsda.srbct <- splsda(X=gooddata.log[,-1], Y=as.factor(gooddata.log$Group), ncomp = 2, keepX = c(100, 100))
plotIndiv(splsda.srbct,ind.names=as.factor(gooddata.log$Group),legend=TRUE)
selectVar(splsda.srbct)

