library(impute)
library(limma)
library(sva)
library(dplyr)


adrc_pos<-read.csv("adrc_pos.csv",check.names=FALSE)
adrc_neg<-read.csv("adrc_neg.csv",check.names=FALSE)
all.equal(adrc_pos$SampleID,adrc_neg$SampleID)

####Read in abundance data and subset to final table of features 
adrc_pos<-read.csv("adrc_pos.csv",check.names=FALSE)[-c(1:2),-c(1:3)]
adrc_pos_mz<-read.csv("adrc_pos.csv",check.names=FALSE)[2,-c(1:3)]
adrc_pos_rt<-read.csv("adrc_pos.csv",check.names=FALSE)[1,-c(1:3)]
adrc_neg<-read.csv("adrc_neg.csv",check.names=FALSE)[-c(1:2),-c(1:3)]
adrc_neg_mz<-read.csv("adrc_neg.csv",check.names=FALSE)[2,-c(1:3)]
adrc_neg_rt<-read.csv("adrc_neg.csv",check.names=FALSE)[1,-c(1:3)]

meta_pheno<-read.csv("adrc_neg.csv",check.names=FALSE)[-c(1:2),c(1:3)]
meta_abunds<-cbind(rbind(adrc_pos_mz,adrc_pos_rt),rbind(adrc_neg_mz,adrc_neg_rt))
abunds<-cbind(adrc_pos,adrc_neg)
mode<-as.factor(c(rep("positive",dim(adrc_pos)[2]),rep("negative",dim(adrc_neg)[2])))
abunds[abunds==0]<-NA
index<-caret::nearZeroVar(abunds)
meta_abunds<-meta_abunds[,-index]
Mode<-mode[-index]
abunds<-abunds[,-index]

####Impute missing data with KNN imputation
library(impute)
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
edata<-log2(apply(edata,2,as.numeric))
rownames(edata)<-rownames(T_abunds)

feature_meta<-cbind(as.data.frame(colnames(abunds)),as.data.frame(Mode),as.data.frame(t(meta_abunds)))
feature_meta_use<-cbind(rownames(T_abunds),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

####Conduct SVA over samples ****method=LEEK RETURNS 0 SVs!!!!!!!!****
pheno_table<-meta_pheno[,1:3]
colnames(pheno_table)<-c("StudyID","DiseaseState","ADRC_CODE")
mod<-model.matrix(~as.factor(DiseaseState),data=pheno_table)
mod0<-model.matrix(~1,data=pheno_table)
n.sv_BE<-num.sv(edata,mod,seed=122)
n.sv_LEEK<-num.sv(edata,mod,seed=122,method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv
total_LEEK<-cbind(as.data.frame(pheno_table),as.data.frame(svobj_LEEK))
total_BE<-cbind(as.data.frame(pheno_table),as.data.frame(svobj_BE))
colnames(total_LEEK)[1:2]<-c("SampleID","Main")
colnames(total_BE)[1:2]<-c("SampleID","Main")
design_table_BE<-cbind(model.matrix(~0+as.factor(Main),data=total_BE),as.data.frame(svobj_BE))
colnames(design_table_BE)[1:4]<-c("Converter","Control")

####Fit linear model for pairwise contrasts 
arrayw<-arrayWeights(edata, design=design_table_BE)  ########design_table_LEEK
fit1<-lmFit(edata,design_table_BE,weights=arrayw)
cm1 <- makeContrasts(`Converter-Control`= Converter-Control,
					 levels=design_table_BE)

fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)


####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="positive")%>%select("MZ","RT","P.Value","Converter.Control"),
			sep="\t",file="ADRC_POS_10_30_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="negative")%>%select("MZ","RT","P.Value","Converter.Control"),
			sep="\t",file="ADRC_NEG_10_30_BE.txt",row.names=FALSE)

cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/ADRC 
mummichog -f ADRC_POS_10_30_BE.txt -o ADRC_POS_10_30_BE -m positive 
mummichog -f ADRC_NEG_10_30_BE.txt -o ADRC_NEG_10_30_BE -m negative 

#########PIUMet
write.table(joinT%>%mutate(`Prize`= -log10(p.adjust(P.Value,method="fdr")))%>%select("MZ","Mode","Prize"),
			sep="\t",file="ADRC_PIUmet_5_2022.txt",row.names=FALSE)


