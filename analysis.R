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
library(Hmisc)

# read and transpose data
#mydata <- read_excel(path="H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\raw data for metaboanalyst no null.xlsx",sheet=1 )
# read in new dataset from Haseeb with identified unknowns
mydata <- read_excel(path="H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\raw data for metaboanalyst no null.xlsx",sheet=1 )

mydata <- as.data.frame(mydata)
mydata <- mydata[!is.na(mydata$Name),]
alldata <- as.data.frame(t(mydata[,-1]))
colnames(alldata) <- mydata$Name

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
gooddata.format_nounk <- gooddata.format[, -grep("UNK",colnames(gooddata.format))]

# impute missing data
nomiss <- MissingValues(gooddata.format,column.cutoff = 0.95,group.cutoff = 0.7,saveoutput = TRUE,
                        outputname = "C:\\Temp\\newoutput",complete.matrix = TRUE)
nomissdf <- fread("C:\\Temp\\newoutput.csv",header=TRUE)
nomissdf <- nomissdf[,-1]
for (i in 1:nrow(nomissdf)) {
  row.names(nomissdf)[i] <- row.names(gooddata.format)[i] 
}
nomissdf <- as.data.frame(nomissdf)
nomissdf_nounk <-  nomissdf[, -grep("UNK",colnames(nomissdf))]
nomissdf.log <- LogTransform(nomissdf)$output

# log transform gooddata
gooddata.log <- LogTransform(gooddata.format)$output
gooddata.log_nounk <- gooddata.log[, -grep("UNK",colnames(gooddata.log))]

#Separating by diet/no diet
gooddata.group<-factor(gooddata.log[,1],levels=unique(gooddata.log[,1]))
dietmat<-gooddata.log[which(gooddata.log[,1]==1),-1]
nodietmat<-gooddata.log[which(gooddata.log[,1]==0),-1]

# creat dataframes with unknowns removed
dietmat_nounk <- dietmat[, -grep("UNK",colnames(dietmat))]
nodietmat_nounk <- nodietmat[, -grep("UNK",colnames(nodietmat))]

#Linear model fit with ordinary statistics
ordFit<-LinearModelFit(datamat=data.matrix(dietmat-nodietmat),
                       ruv2=FALSE,
                       factormat=matrix(1,nrow=nrow(dietmat)),
                       outputname = "H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\ordFit",
                       saveoutput = TRUE)
TwoGroupPlots(gooddata.log[,-1],
              tstats = ordFit$t[,1],
              foldchanges = ordFit$coef[,1],
              pvalues = ordFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)
# some compounds have NA's for p-values

#Linear model fit with moderated statistics
modFit<-LinearModelFit(datamat=data.matrix(dietmat-nodietmat),
                       ruv2=FALSE,
                       moderated=TRUE,
                       factormat=matrix(1,nrow=nrow(dietmat)),
                       outputname = "H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\modFit",
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

#Linear model fit with moderated statistics - exclude unknowns
modFit<-LinearModelFit(datamat=data.matrix(dietmat_nounk-nodietmat_nounk),
                       ruv2=FALSE,
                       moderated=TRUE,
                       factormat=matrix(1,nrow=nrow(dietmat)),
                       outputname = "H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\modFit_nounk",
                       saveoutput = TRUE)
# get fold change for knowns
fc_nounk <- as.data.frame(FoldChange(gooddata.log_nounk,paired=TRUE))
fc_nounk$Compound <- rownames(fc_nounk)
fc_nounk <- as.data.frame(fc_nounk[,-1])
colnames(fc_nounk) <- c("FC","X")
# merge fold change with results of moderated t-test for input into Metscape
temp <- read.csv("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\modFit_nounk.csv")
metscape <- merge(temp,fc_nounk,by="X")
write.csv(metscape,"H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\formetscape.csv")

# Volcano plot
modfit_nona <- modFit[!is.na(modFit$p.value),]
VolcanoPlot(folds=modfit_nona$coef, pvals=modfit_nona$p.value,plimit=0.05)
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

# looking at direction of differences
MetBoxPlots(gooddata.log,"ISOLEUCINE (M+H)+")
MetBoxPlots(gooddata.log,"ARACHIDONIC ACID (M-H)-")

# Dendrogram
Dendrogram(gooddata.log_nounk)
# is there an interaction between diet and PCOS?  PCO group seem like they are more similar regardless of diet

# HeatMap 
HeatMap(nomissdf_nounk,colramp=redgreen(75),margins = c(5,10),key=FALSE,dendrogram = "both")

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
probject <- pcaMethods::pca(md,method="nipals",nPcs = 3)
plotPcs(probject,pcs=1:3,type="scores",col=as.factor(gooddata.format$Group))

# create dataset for PLS-DA
# create variable for PCO status
gooddata.plsda <- gooddata.format
for (i in 1:nrow(gooddata.plsda)) {
  gooddata.plsda$PCOS[i] <- ifelse(substring(row.names(gooddata.plsda[i,]),1,1)=="P",1,0)
}
# reorder dataset
gooddata.plsda <- gooddata.plsda[,c(1,6690,2:6689)]
# do the same for the nomiss dataset
nomissdf <- fread("C:\\Temp\\newoutput.csv",header=TRUE)
nomissdf <- nomissdf[,-1]
for (i in 1:nrow(nomissdf)) {
  row.names(nomissdf)[i] <- row.names(gooddata.format)[i] 
}
nomiss.plsda <- nomissdf
nomiss.plsda$PCOS <- substring(row.names(nomiss.plsda),1,1)
nomiss.plsda$id <- row.names(nomiss.plsda)
nomiss.plsda$id <- gsub("P", "", nomiss.plsda$id)
nomiss.plsda$id <- gsub("C", "", nomiss.plsda$id)
nomiss.plsda$id <- gsub("diet", "", nomiss.plsda$id)
nomiss.plsda$id <- gsub("No-", "", nomiss.plsda$id)

# create dataset of demographics
demo <- gooddata.plsda[c("PCOS")]
demo$id <- row.names(nomiss.plsda)
demo$id <- gsub("P", "", nomiss.plsda$id)
demo$id <- gsub("C", "", nomiss.plsda$id)
demo$id <- gsub("diet", "", nomiss.plsda$id)
demo$id <- gsub("No-", "", nomiss.plsda$id)
temp <- read.csv("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\demographics.csv")
temp <- temp[c("subject_id","age","gender","ethnicity","tanner","bmi_percentile")]
temp$id <- gsub("6164-", "", temp$subject_id) 
demo <- merge(demo,temp,by="id")
demo <- demo %>% distinct(id, .keep_all=TRUE)
demo$dummy <- rep(1,nrow(demo))
demo$gender <- as.factor(demo$gender)
demo$tanner <- as.factor(demo$tanner)
demo$bmi_percentile <- as.numeric(as.character(demo$bmi_percentile))
label(demo$age)="Age"
label(demo$gender)="Gender"
label(demo$ethnicity)="Ethnicity"
label(demo$tanner)="Tanner"
label(demo$bmi_percentile)="BMI %ile"

# table 1
tab1 <- final_table(data=demo,variables=c("age","gender","ethnicity","tanner","bmi_percentile"),
                    ron=2,group=as.factor(demo$dummy),margin=2)

# http://mixomics.org/mixmc/case-study-hmp-bodysites-repeated-measures/
splsda.diet = splsda(X = nomiss.plsda[,-c(1,6689:6690)], Y=as.factor(nomiss.plsda$Group), 
                   ncomp = 2, multilevel = as.factor(nomiss.plsda$id),keepX = c(200, 200))
plotIndiv(splsda.diet, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = FALSE, title="")
listvar <- selectVar(splsda.diet)
listvar <- listvar$name[-grep("UNK",listvar$name)]
set.seed(34)  # for reproducible results for this code
diet.perf.splsda = perf(splsda.diet, validation = 'Mfold', folds = 5, 
                           progressBar = FALSE, nrepeat = 10, dist = 'max.dist',auc=TRUE)
diet.perf.splsda$error.rate
plot(diet.perf.splsda)
head(selectVar(splsda.diet, comp = 1)$value) 
cim(splsda.diet, row.sideColors = color.mixo(as.factor(nomiss.plsda$Group)))
diet.perf.splsda.loo = perf(splsda.diet, validation = 'loo', 
                        progressBar = FALSE, auc=TRUE)
diet.auroc <- auroc(splsda.diet)

# biplot with top 20 compounds
splsda.diet20 = splsda(X = nomiss.plsda[,-c(1,6689:6690)], Y=as.factor(nomiss.plsda$Group), 
                     ncomp = 2, multilevel = as.factor(nomiss.plsda$id),keepX = c(20, 20))
ind.coord <- splsda.diet20$variates$X[, 1:2]
var.coord = plotVar(splsda.diet20,var.names = FALSE)[,c("x","y")]
biplot(ind.coord,var.coord,xlabs=as.factor(nomiss.plsda$Group))
abline(h=0,v=0,lty=2)

# boxplots of known compounds from the top 200 from PLSDA
MetBoxPlots(gooddata.log,"2,3-BISPHOSPHO-D-GLYCERIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"3,4-DIHYDROXYBENZOIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"TAURINE (M-H)-")
MetBoxPlots(gooddata.log,"INOSINE 5'-DIPHOSPHATE (M-H)-")
MetBoxPlots(gooddata.log,"MEVALONIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"BENZOIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"DETHIOBIOTIN (M+H)+[-H2O]")
MetBoxPlots(gooddata.log,"O-SUCCINYL-L-HOMOSERINE (M-H)-")
MetBoxPlots(gooddata.log,"CYSTEIC ACID (M-H)-")
MetBoxPlots(gooddata.log,"DIHYDROXYFUMARIC ACID (M-H)-")

# now for pcos
# will not converge
splsda.pcos = splsda(X = nomiss.plsda[,-c(1,6689:6690)], Y=as.factor(nomiss.plsda$PCOS), 
                     ncomp = 2, multilevel = as.factor(nomiss.plsda$id),max.iter = 10000,
                     keepX = c(10, 10))
a <- tune.splsda(X = nomiss.plsda[,-c(1,6689:6690)], Y=as.factor(nomiss.plsda$PCOS), 
            ncomp = 2, multilevel = as.factor(nomiss.plsda$id))
plotIndiv(splsda.pcos, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE)
set.seed(34)  # for reproducible results for this code
pcos.perf.splsda = perf(splsda.pcos, validation = 'Mfold', folds = 5, 
                        progressBar = FALSE, nrepeat = 10, dist = 'max.dist',auc=TRUE)
pcos.perf.splsda$error.rate
plot(pcos.perf.splsda)
head(selectVar(splsda.pcos, comp = 1)$value) 
cim(splsda.pcos, row.sideColors = color.mixo(as.factor(nomiss.plsda$PCOS)))
pcos.auroc <- auroc(splsda.pcos)

# ROC for splsda
auroc(splsda.srbct,newdata=splsda.srbct$input.X,outcome.test = as.factor(splsda.srbct$Y),plot=TRUE)
auroc(splsda.pcos,newdata=splsda.pcos$input.X,outcome.test = as.factor(splsda.pcos$Y),plot=TRUE)

