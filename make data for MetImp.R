
setwd("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\")

library(readxl)
# load NormalizeMets first, metabolomics second
library(NormalizeMets)
#library(metabolomics)
library(pca3d)
library(pcaMethods)
library(ggplot2)
library(rgl)
library(stringr)
library(data.table)
library(dplyr)
library(mixOmics)
library(Hmisc)
library("FactoMineR")
library("factoextra")
library(tableone)
library(gdata)
library(forcats)
library(limma)

# read and transpose data
#mydata <- read_excel(path="H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\raw data for metaboanalyst no null.xlsx",sheet=1 )
# read in new dataset from Haseeb with identified unknowns
mydata <- read_excel(path="H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\updated data with identified unknowns 1.2.19 no null.xlsx",sheet=1 )
mydata <- as.data.frame(mydata)
mydata <- mydata[!is.na(mydata$Name),]
alldata <- as.data.frame(t(mydata[,-1]))
colnames(alldata) <- mydata$Name
colnames(alldata) <- str_remove_all(colnames(alldata),fixed(" (M-H)-"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+Cl)-"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M-H)-[-H2O]"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M-H)+[-H2O]"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+Cl)+[-H2O]"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+Cl)-"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+H)-"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" [-H2O]"))
colnames(alldata) <- str_remove(colnames(alldata),fixed("[-H2O]"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (2M-H)+"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+H)+"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+K)+"))
colnames(alldata) <- str_remove(colnames(alldata),fixed(" (M+Na)+"))

source("C:\\Users\\pylell\\Documents\\GitHub\\General-code\\temp_table1.r")
source("C:\\Users\\pylell\\Documents\\GitHub\\General-code\\foldchange.r")
source("C:\\Users\\pylell\\Documents\\GitHub\\General-code\\editcolnames.r")

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

# check how  many compounds are known
check_knowns <- alldata[, -grep("UNK",colnames(alldata))]

# convert to numeric
num <- alldata[ ,!(colnames(alldata) %in% c("Batch","id","longid","uniqueid"))]
fwrite(num,"C:\\Temp\\temp.csv")
newnum <- fread("C:\\Temp\\temp.csv",colClasses = "numeric")
alldata <- alldata[ ,(colnames(alldata) %in% c("Batch","id","uniqueid"))]
alldata <- cbind(alldata,newnum)

# keep only diet condition since that has less missing data
gooddata <- alldata[alldata$Batch=="diet",]
for (i in 1:nrow(gooddata)) {
  gooddata$PCOS[i] <- ifelse(substring(row.names(gooddata[i,]),1,1)=="P",1,0)
}
gooddata$group <- gooddata$PCOS
gooddata <- cbind(gooddata[,1:3],gooddata[,6748:6749],gooddata[,4:6747])
gooddata$Batch <- NULL
gooddata$id <- NULL
gooddata$PCOS <- NULL

write.csv(gooddata,"H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\for_metimp.csv")