library(dplyr)
library(mixOmics)

G<-read.csv("Complete_Pheno.csv")
R<-read.csv("join_quanterix.csv")
colnames(R)[1]<-"OtherID"
U<-dplyr::inner_join(R,G,by="OtherID")
U<-U%>%mutate(`AB42/AB40`=U$`Conc_AB42`/U$`Conc_AB40`)
U<-U%>%mutate(`AB40/AB42`=U$`Conc_AB40`/U$`Conc_AB42`)
J<-pca(U[,c(2:5,29,30)],ncomp=2)
plotIndiv(J,group=U$Synd,legend=TRUE)
quantile(J$variates$X[,1],probs=seq(0,1,.1))

wilcox.test(J$variates$X[,1]~as.factor(U$Synd))
vioplot(J$variates$X[,1]~as.factor(U$Synd))

k<-U[,c(2:5,29,30)]
vioplot(k[,]~as.factor(U$Synd))
wilcox.test(k[,]~as.factor(U$Synd))


#############################
#############################
G<-read.csv("Complete_Pheno.csv")
R<-read.csv("join_Quanterix_NfL_2More.csv")
colnames(R)[1]<-"OtherID"
U<-dplyr::inner_join(R,G,by="OtherID")
wilcox.test(U[,5]~as.factor(U$Synd))
vioplot(U[,5]~as.factor(U$Synd))




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

pheno_table<-meta_pheno[,1:3]
colnames(pheno_table)<-c("StudyID","DiseaseState","ADRC_CODE")











G<-read.csv("Complete_Pheno.csv")
R<-read.csv("join_quanterix.csv")
colnames(R)[1]<-"OtherID"
U<-dplyr::inner_join(R,G,by="OtherID")
corData<-cbind(pheno_table,as.data.frame(t(edata)))
colnames(corData)[1]<-"OtherID"


joinCor<-dplyr::inner_join(U,corData,by="OtherID")
reduced_Preclin<-joinCor[,c(1:3,10,31:dim(joinCor)[2])]%>%filter(Synd=="Converterpre")
reduced_Control<-joinCor[,c(1:3,10,31:dim(joinCor)[2])]%>%filter(Synd=="normal")

ab40_cor_Preclin<- numeric(ncol(reduced_Preclin[,-c(1:4)])) 
ab40_P.Value_Preclin<- numeric(ncol(reduced_Preclin[,-c(1:4)]))
for(i in 5:ncol(reduced_Preclin[,-c(1:4)])){
         ab40_cor_Preclin[i]<-cor.test(reduced_Preclin[,i],reduced_Preclin$Conc_AB40,method="spearman",exact=FALSE)$estimate
         ab40_P.Value_Preclin[i]<-cor.test(reduced_Preclin[,i],reduced_Preclin$Conc_AB40,method="spearman",exact=FALSE)$p.value
}
A<-cbind(feature_meta_use,ab40_cor_Preclin,ab40_P.Value_Preclin)
colnames(A)[c(4,5,7,6)]<-c("MZ","RT","P.Value","ab40_cor_Preclin")

write.table(A%>%filter(Mode=="positive")%>%dplyr::select("MZ","RT","P.Value","ab40_cor_Preclin"),
			sep="\t",file="ADRC_CORR_4_7_AB40_Preclin.POS_MUMMI.txt",row.names=FALSE)
write.table(A%>%filter(Mode=="negative")%>%dplyr::select("MZ","RT","P.Value","ab40_cor_Preclin"),
			sep="\t",file="ADRC_CORR_4_7_AB40_Preclin.NEG_MUMMI.txt",row.names=FALSE)

cd /Volumes/TJGross_Remote/318557/Active_Projects_4_9_2021/ADRC

mummichog -f ADRC_CORR_4_7_AB40_Preclin.POS_MUMMI.txt -o ADRC_CORR_4_7_AB40_Preclin.POS_MUMMI -m positive
mummichog -f ADRC_CORR_4_7_AB40_Preclin.NEG_MUMMI.txt -o ADRC_CORR_4_7_AB40_Preclin.NEG_MUMMI -m negative 

write.table(A%>%dplyr::mutate(`Prize`= -log10(P.Value))%>%dplyr::select("MZ","Mode","Prize"),
			sep="\t",file="ADRC_CORR_4_7_PIUmet_AB40_Preclin.txt",row.names=FALSE)


ab40_cor_Control<- numeric(ncol(reduced_Control[,-c(1:4)])) 
ab40_P.Value_Control<- numeric(ncol(reduced_Control[,-c(1:4)]))
for(i in 5:ncol(reduced_Control[,-c(1:4)])){
         ab40_cor_Control[i]<-cor.test(reduced_Control[,i],reduced_Control$Conc_AB40,method="spearman",exact=FALSE)$estimate
         ab40_P.Value_Control[i]<-cor.test(reduced_Control[,i],reduced_Control$Conc_AB40,method="spearman",exact=FALSE)$p.value
}
B<-cbind(feature_meta_use,ab40_cor_Control,ab40_P.Value_Control)
colnames(B)[c(4,5,7,6)]<-c("MZ","RT","P.Value","ab40_cor_Control")

write.table(B%>%filter(Mode=="positive")%>%dplyr::select("MZ","RT","P.Value","ab40_cor_Control"),
			sep="\t",file="ADRC_CORR_4_7_AB40_Control.POS_MUMMI.txt",row.names=FALSE)
write.table(B%>%filter(Mode=="negative")%>%dplyr::select("MZ","RT","P.Value","ab40_cor_Control"),
			sep="\t",file="ADRC_CORR_4_7_AB40_Control.NEG_MUMMI.txt",row.names=FALSE)

cd /Volumes/TJGross_Remote/318557/Active_Projects_4_9_2021/ADRC

mummichog -f ADRC_CORR_4_7_AB40_Control.POS_MUMMI.txt -o ADRC_CORR_4_7_AB40_Control.POS_MUMMI -m positive
mummichog -f ADRC_CORR_4_7_AB40_Control.NEG_MUMMI.txt -o ADRC_CORR_4_7_AB40_Control.NEG_MUMMI -m negative 


B[B$P.Value<0.00000000]<-0.000000001
write.table(B%>%dplyr::mutate(`Prize`= -log10(P.Value))%>%dplyr::select("MZ","Mode","Prize"),
			sep="\t",file="ADRC_CORR_4_7_PIUmet_AB40_Control.txt",row.names=FALSE)
