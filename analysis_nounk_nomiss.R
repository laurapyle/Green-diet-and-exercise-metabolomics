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
nomissdf_nounk.log <- LogTransform(nomissdf_nounk)$output

#Separating by diet/no diet
nomissdf_nounk.group<-factor(nomissdf_nounk.log[,1],levels=unique(nomissdf_nounk.log[,1]))
dietmat_nounk<-nomissdf_nounk.log[which(nomissdf_nounk.log[,1]==1),-1]
nodietmat_nounk<-nomissdf_nounk.log[which(nomissdf_nounk.log[,1]==0),-1]

#Linear model fit with ordinary statistics
ordFit<-LinearModelFit(datamat=data.matrix(dietmat_nounk-nodietmat_nounk),
                       ruv2=FALSE,
                       factormat=matrix(1,nrow=nrow(dietmat_nounk)),
                       outputname = "H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\ordFit",
                       saveoutput = TRUE)
TwoGroupPlots(nomissdf_nounk.log[,-1],
              tstats = ordFit$t[,1],
              foldchanges = ordFit$coef[,1],
              pvalues = ordFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)
# some compounds have NA's for p-values

#Linear model fit with moderated statistics
modFit<-LinearModelFit(datamat=data.matrix(dietmat_nounk-nodietmat_nounk),
                       ruv2=FALSE,
                       moderated=TRUE,
                       factormat=matrix(1,nrow=nrow(dietmat_nounk)),
                       outputname = "H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\modFit",
                       saveoutput = TRUE)
TwoGroupPlots(nomissdf_nounk.log[,-1],
              tstats = modFit$t[,1],
              foldchanges = modFit$coef[,1],
              pvalues = modFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(2),
              pcutoff = 0.05)
# trying TwoGroupPlot with anything significantly different, regardless of fold change
TwoGroupPlots(nomissdf_nounk.log[,-1],
              tstats = modFit$t[,1],
              foldchanges = modFit$coef[,1],
              pvalues = modFit$p.val[,1],
              padjmethod = "BH",
              fcutoff = log(1),
              pcutoff = 0.05)

# get fold change for knowns
fc_nounk <- as.data.frame(FoldChange(nomissdf_nounk.log,paired=TRUE))
fc_nounk$Compound <- rownames(fc_nounk)
fc_nounk <- as.data.frame(fc_nounk[,-1])
colnames(fc_nounk) <- c("FC","X")
# merge fold change with results of moderated t-test for input into Metscape
temp <- read.csv("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\modFit_nounk.csv")
metscape <- merge(temp,fc_nounk,by="X")
metscape <- metscape[,c(1,8,10)]
annotated <- read.csv("H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\formetscape_clean.csv")
annotated <- annotated[,c(1:2)]
metscape <- merge(metscape,annotated,by="X")
write.csv(metscape,"H:\\Endocrinology\\Green\\Metabolomics papers\\Diet and exercise\\Data\\formetscape_imputed.csv")

# Volcano plot
modfit_nona <- modFit[!is.na(modFit$p.value),]
VolcanoPlot(folds=modfit_nona$coef, pvals=modfit_nona$p.value,plimit=0.05)
# VolcanoPlot doesn't like NA's for p-values

# Dendrogram
Dendrogram(nomissdf_nounk.log)
# is there an interaction between diet and PCOS?  PCO group seem like they are more similar regardless of diet

# HeatMap 
HeatMap(nomissdf_nounk.log,colramp=redgreen(75),margins = c(5,10),key=FALSE,dendrogram = "both")

# PCA plot
PcaPlots(nomissdf_nounk.log,scale=TRUE, center=TRUE)

# RLA plot
RlaPlots(nomissdf_nounk.log,type="ag")
RlaPlots(nomissdf_nounk.log,type="wg")

# create dataset for PLS-DA
# create variable for PCO status
nomiss.plsda <- nomissdf_nounk
for (i in 1:nrow(nomissdf_nounk)) {
  nomiss.plsda$PCOS[i] <- ifelse(substring(row.names(nomissdf_nounk[i,]),1,1)=="P",1,0)
}
nomiss.plsda$id <- row.names(nomiss.plsda)
nomiss.plsda$id <- gsub("P", "", nomiss.plsda$id)
nomiss.plsda$id <- gsub("C", "", nomiss.plsda$id)
nomiss.plsda$id <- gsub("diet", "", nomiss.plsda$id)
nomiss.plsda$id <- gsub("No-", "", nomiss.plsda$id)
nomiss.plsda <- nomiss.plsda[,c(1,167,168,2:166)]

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
splsda.diet = splsda(X = nomiss.plsda[,c(4:168)], Y=as.factor(nomiss.plsda$Group), 
                   ncomp = 2, multilevel = as.factor(nomiss.plsda$id))
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
plot.new()
cim(splsda.diet, row.sideColors = color.mixo(as.factor(nomiss.plsda$Group)))
diet.perf.splsda.loo = perf(splsda.diet, validation = 'loo', 
                        progressBar = FALSE, auc=TRUE)
diet.auroc <- auroc(splsda.diet)

# biplot with top 20 compounds
splsda.diet20 = splsda(X = nomiss.plsda[,c(4:168)], Y=as.factor(nomiss.plsda$Group), 
                     ncomp = 2, multilevel = as.factor(nomiss.plsda$id),keepX = c(20, 20))
ind.coord <- splsda.diet20$variates$X[, 1:2]
var.coord = plotVar(splsda.diet20,var.names = FALSE)[,c("x","y")]
biplot(ind.coord,var.coord,xlabs=as.factor(nomiss.plsda$Group))
abline(h=0,v=0,lty=2)

# now for pcos
# need to do imputation by PCOS group
gooddata.pcos <- gooddata.format
for (i in 1:nrow(gooddata.pcos)) {
  gooddata.pcos$PCOS[i] <- ifelse(substring(row.names(gooddata.pcos[i,]),1,1)=="P",1,0)
}
gooddata.pcos$Group <- gooddata.pcos$PCOS
nomiss.pcos <- MissingValues(gooddata.pcos,column.cutoff = 0.95,group.cutoff = 0.7,saveoutput = TRUE,
                        outputname = "C:\\Temp\\newoutput_pcos",complete.matrix = TRUE)
nomissdf_pcos <- fread("C:\\Temp\\newoutput_pcos.csv",header=TRUE)
nomissdf_pcos <- nomissdf_pcos[,-1]
for (i in 1:nrow(nomissdf_pcos)) {
  row.names(nomissdf_pcos)[i] <- row.names(gooddata.pcos)[i] 
}
nomissdf_pcos <- as.data.frame(nomissdf_pcos)
nomissdf_nounk_pcos <-  nomissdf[, -grep("UNK",colnames(nomissdf_pcos))]
for (i in 1:nrow(nomissdf_nounk_pcos)) {
  nomissdf_nounk_pcos$PCOS[i] <- ifelse(substring(row.names(nomissdf_nounk_pcos[i,]),1,1)=="P",1,0)
}
nomissdf_nounk_pcos$Group <- nomissdf_nounk_pcos$PCOS
nomissdf_nounk_pcos.log <- LogTransform(nomissdf_nounk_pcos)$output
nomiss.plsda_pcos <- nomissdf_nounk_pcos
nomiss.plsda_pcos$id <- row.names(nomiss.plsda_pcos)
nomiss.plsda_pcos$id <- gsub("P", "", nomiss.plsda_pcos$id)
nomiss.plsda_pcos$id <- gsub("C", "", nomiss.plsda_pcos$id)
nomiss.plsda_pcos$id <- gsub("diet", "", nomiss.plsda_pcos$id)
nomiss.plsda_pcos$id <- gsub("No-", "", nomiss.plsda_pcos$id)
nomiss.plsda_pcos <- nomiss.plsda_pcos[,c(1,167,168,2:166)]
splsda.pcos = splsda(X = nomiss.plsda_pcos[,c(4:168)], Y=as.factor(nomiss.plsda_pcos$PCOS), 
                     ncomp = 2, multilevel = as.factor(nomiss.plsda_pcos$id),max.iter = 10000)
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
# why is this not converging?
# probably because there are only a handful of controls
nomiss.plsda_pcos.log <- nomiss.plsda_pcos[,-(2:3)]
nomiss.plsda_pcos.log <- LogTransform(nomiss.plsda_pcos.log)$output
# HeatMap 
HeatMap(nomiss.plsda_pcos.log,colramp=redgreen(75),margins = c(5,10),key=FALSE,dendrogram = "both")
# PCA plot
PcaPlots(nomiss.plsda_pcos.log,scale=TRUE, center=TRUE)

# can we do a PCA plot with 4 groups (i.e., cross-tab of PCO and diet)
fourgroup <- nomiss.plsda
fourgroup$temp[fourgroup$PCOS==0 & fourgroup$Group==0] <- 1
fourgroup$temp[fourgroup$PCOS==0 & fourgroup$Group==1] <- 2
fourgroup$temp[fourgroup$PCOS==1 & fourgroup$Group==0] <- 3
fourgroup$temp[fourgroup$PCOS==1 & fourgroup$Group==1] <- 4
fourgroup <- fourgroup[,c(1:3,169,4:168)]
fourgroup$Group <- fourgroup$temp
fourgroup <- fourgroup[,-(2:4)]
fourgroup.log <- LogTransform(fourgroup)$output
# HeatMap 
HeatMap(fourgroup.log,colramp=redgreen(75),margins = c(5,10),key=FALSE,dendrogram = "both")
# PCA plot
PcaPlots(fourgroup.log,scale=TRUE, center=TRUE)
# this is really cool - the diet and nodiet cluster without respect to PCOS status

