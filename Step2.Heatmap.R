rm(list=ls())
library(dplyr)
library(gplots)
library(survival)
library(rpart)


high_grade <- read.csv("C:/Users/korea/Desktop/COAD/clinical/high_grade_colon_mRNA.csv", row.names = 1)
low_grade <- read.csv("C:/Users/korea/Desktop/COAD/clinical/low_grade_colon_mRNA.csv", row.names = 1)

high_grade <- log2(high_grade+1)
low_grade <- log2(low_grade+1)

#######heatmap - wilcox
selected_gene_high <- high_grade[,c("TEAD3","C6orf15","ARHGDIG","B3GALT4","ITFG3","FKBP9","MC1R","BAT3","TCHH","MYST1",
                                    "RARA","CAMK2B","PCDHA2","CARKD","DDAH2","NPR3","PHLDB3","TCF7L1","KIAA0174","WIPI2",
                                    "RGL2","FBRS","HAP1","TMEM198","CDIPT","CAPS","C13orf15","ZNF821","ATF6B","SHC1",
                                    "FADS6","EEF1A2","LRCH4","STXBP5L","LYPD3","C6orf106","AGPAT1","PPL","SCN5A","LAMP1",
                                    "TMEM105","C5orf23","RIMS4","GLYR1","LEMD2","PCGF2","TERF2IP","UPK2","LRP2","PALM2",
                                    "AGPAT5","COQ2","PDE12","NAT1","PBK","ERI1","INTS10","GSR","NUDT6","LYSMD2",
                                    "MINPP1","CDCA2","BEND3","ME2","MCPH1","ESCO2","CHAC2","MDM2","PINX1","CASP1",
                                    "CXCL11","TNFRSF11A","C12orf11","CNOT7","UNG","FAS","NEIL2","NEK4","APOL6","FASLG",
                                    "ACADSB","GTF2E2","FUT10","XPO7","FAM83F","C18orf55","MMAA","ATP6V1B2","C2CD4A","EXD2",
                                    "TNFRSF10A","COX10","YARS2","C12orf48","DDX20","PPA2","AMD1","SCO1","STX18","MORC4"
)]
selected_gene_low <- low_grade[,c("TEAD3","C6orf15","ARHGDIG","B3GALT4","ITFG3","FKBP9","MC1R","BAT3","TCHH","MYST1",
                                  "RARA","CAMK2B","PCDHA2","CARKD","DDAH2","NPR3","PHLDB3","TCF7L1","KIAA0174","WIPI2",
                                  "RGL2","FBRS","HAP1","TMEM198","CDIPT","CAPS","C13orf15","ZNF821","ATF6B","SHC1",
                                  "FADS6","EEF1A2","LRCH4","STXBP5L","LYPD3","C6orf106","AGPAT1","PPL","SCN5A","LAMP1",
                                  "TMEM105","C5orf23","RIMS4","GLYR1","LEMD2","PCGF2","TERF2IP","UPK2","LRP2","PALM2",
                                  "AGPAT5","COQ2","PDE12","NAT1","PBK","ERI1","INTS10","GSR","NUDT6","LYSMD2",
                                  "MINPP1","CDCA2","BEND3","ME2","MCPH1","ESCO2","CHAC2","MDM2","PINX1","CASP1",
                                  "CXCL11","TNFRSF11A","C12orf11","CNOT7","UNG","FAS","NEIL2","NEK4","APOL6","FASLG",
                                  "ACADSB","GTF2E2","FUT10","XPO7","FAM83F","C18orf55","MMAA","ATP6V1B2","C2CD4A","EXD2",
                                  "TNFRSF10A","COX10","YARS2","C12orf48","DDX20","PPA2","AMD1","SCO1","STX18","MORC4"
)]
total_gene_expression <- rbind(selected_gene_low, selected_gene_high)

gene_values<-total_gene_expression
Samples_mean<-apply(gene_values,2,mean)   # Standardize the data
Samples_sd<-apply(gene_values,2,sd)
Samples_stand<-gene_values
for(i in 1:dim(gene_values)[2])
{
  Samples_stand[,i]<-(gene_values[,i]-Samples_mean[i])/Samples_sd[i]
}

total_gene_expression<-Samples_stand

select_gene_low <- total_gene_expression[c(1:nrow(selected_gene_low)), ]
select_gene_high <- total_gene_expression[c(nrow(selected_gene_low)+1:nrow(selected_gene_high)), ]

select_gene_low_trans <- t(select_gene_low)
select_gene_high_trans <- t(select_gene_high)

select_gene_low_trans_arrange <- select_gene_low_trans[order(apply(select_gene_low_trans, 1, mean) - apply(select_gene_high_trans, 1, mean)), ]
select_gene_high_trans_arrange <- select_gene_high_trans[order(apply(select_gene_low_trans, 1, mean) - apply(select_gene_high_trans, 1, mean)), ]

combinded_heatmap_data<-cbind(select_gene_low_trans,NA,NA,NA, select_gene_high_trans)

combinded_heatmap_data_2<-matrix(0,nrow=dim(combinded_heatmap_data)[1],ncol=dim(combinded_heatmap_data)[2])
for(i in 1:dim(combinded_heatmap_data)[2])
{
  combinded_heatmap_data_2[,i]<-as.numeric(combinded_heatmap_data[,i])
}

rownames(combinded_heatmap_data_2)<-rownames(combinded_heatmap_data)
col_names<-colnames(combinded_heatmap_data)
col_names[col_names=="NA"]<-""
colnames(combinded_heatmap_data_2)<-col_names

file_name_plots=paste("C:/Users/korea/Desktop/COAD/step2.heatmap/heatmap_0407.eps", sep="")
postscript(file=file_name_plots, paper="a4",horizontal =TRUE) 
par(mfrow=c(2,1) )
par(mar=c(2.1,1.1,2.1,1.1))
par(oma=c(1,1,1,1))

my_palette <- colorRampPalette(c("green", "black", "red"))(n = 100)
heatmap.2(combinded_heatmap_data_2,col=my_palette,Rowv=FALSE,Colv=FALSE,density.info="none",trace="none", dendrogram="none"
          ,main="Coloncancer_DEG")

dev.off()

