library(sva)
library(dplyr)
adrc_pos<-read.csv("adrc_pos.csv",check.names=FALSE)
adrc_neg<-read.csv("adrc_neg.csv",check.names=FALSE)
identical(adrc_pos$SampleID,adrc_neg$SampleID)

####Read in abundance data and subset to final table of features 
adrc_pos<-read.csv("adrc_pos.csv",check.names=FALSE)[-c(1:2),-c(1:3)]
adrc_pos_mz<-read.csv("adrc_pos.csv",check.names=FALSE)[2,-c(1:3)]
adrc_pos_rt<-read.csv("adrc_pos.csv",check.names=FALSE)[1,-c(1:3)]
adrc_neg<-read.csv("adrc_neg.csv",check.names=FALSE)[-c(1:2),-c(1:3)]
adrc_neg_mz<-read.csv("adrc_neg.csv",check.names=FALSE)[2,-c(1:3)]
adrc_neg_rt<-read.csv("adrc_neg.csv",check.names=FALSE)[1,-c(1:3)]

meta_pheno<-read.csv("adrc_neg.csv",check.names=FALSE)[-c(1:2),c(1:3)]
meta_abunds<-cbind(rbind(adrc_pos_mz,adrc_pos_rt),rbind(adrc_neg_mz,adrc_neg_rt))
abunds<-log(cbind(adrc_pos,adrc_neg))
mode<-as.factor(c(rep("ESI+",dim(adrc_pos)[2]),rep("ESI-",dim(adrc_neg)[2])))
na_index<-apply(abunds,2,is.na)
infinite_index<-apply(abunds,2,is.infinite)
final_index<-na_index==TRUE|infinite_index==TRUE
abunds[final_index==TRUE]<-NA
index<-caret::nearZeroVar(abunds)
meta_abunds<-meta_abunds[,-index]
mode<-mode[-index]
abunds<-abunds[,-index]

####Impute missing data with KNN imputation
library(impute)
T_abunds<-t(abunds)
imputed_log<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)

####Carry out SVA
edata<-imputed_log$data
mod<-model.matrix(~as.factor(factor(pheno)), data=meta_pheno)
mod0<-model.matrix(~1,data=meta_pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)
modSv<-cbind(mod,svobj$sv)
mod0Sv<-cbind(mod0,svobj$sv)
pValuesSv<-f.pvalue(edata,modSv,mod0Sv)
qValuesSv<-p.adjust(pValuesSv,method="fdr")

median_case<-apply(as.data.frame(edata[,meta_pheno$pheno=="Case"]),1,median)
median_control<-apply(edata[,meta_pheno$pheno=="Control"],1,median)
log2FC<-log2(exp(median_case)/exp(median_control))
mz<-gsub("_.+$","",rownames(edata))
rt<-gsub("^.+_","",rownames(edata))
table<-cbind(as.data.frame(t(meta_abunds)),as.data.frame(mode),as.data.frame(pValuesSv),as.data.frame(qValuesSv),as.data.frame(log2FC))
colnames(table)<-c("MZ","RT","Mode","PValue","QValue","Log2FC")
POS_ADRC<-table%>%filter(Mode=="ESI+")
NEG_ADRC<-table%>%filter(Mode=="ESI-")
write.csv(POS_ADRC, file="SVA_DE_ADRC_POS.csv")
write.csv(NEG_ADRC, file="SVA_DE_ADRC_NEG.csv")