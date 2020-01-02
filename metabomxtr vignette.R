library(metabomxtr)
data(metabdata)
yvars<-colnames(metabdata)[24:27]
metabdata$PHENO<-relevel(metabdata$PHENO,ref="MomLowFPG")
fullModel<-~PHENO|PHENO+age_ogtt_mc+ga_ogtt_wks_mc+storageTimesYears_mc+parity12
reducedModel<-~1|age_ogtt_mc+ga_ogtt_wks_mc+storageTimesYears_mc+parity12
fullModelResults<-mxtrmod(ynames=yvars,mxtrModel=fullModel,data=metabdata)
