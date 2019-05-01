rm(list=ls())
library("survival")
library("glmnet")


file_name_tumor=paste("C:/Users/korea/desktop/COAD/clinical/combined_data.csv", sep="")
tumor_Data <- read.table(file_name_tumor, header=TRUE, sep=",")


Day.to.follow.up <- as.numeric(as.character(tumor_Data[,"days_to_last_followup"]))
Day.to.Death <- as.numeric(as.character(tumor_Data[,"days_to_death"]))

##DEG&positive&negative
SLC22A17 <- log2(tumor_Data[,"SLC22A17"]+1)
APBB1 <- log2(tumor_Data[,"APBB1"]+1)
SLC7A14 <- log2(tumor_Data[,"SLC7A14"]+1)
JAM3 <- log2(tumor_Data[,"JAM3"]+1)
PRELP <- log2(tumor_Data[,"PRELP"]+1)
LYSMD3 <- log2(tumor_Data[,"LYSMD3"]+1)
RBPMS2 <- log2(tumor_Data[,"RBPMS2"]+1)
LMOD1 <- log2(tumor_Data[,"LMOD1"]+1)
GEFT <- log2(tumor_Data[,"GEFT"]+1)
SALL2 <- log2(tumor_Data[,"SALL2"]+1)
TNS1 <- log2(tumor_Data[,"TNS1"]+1)
EFEMP2 <- log2(tumor_Data[,"EFEMP2"]+1)
SYDE1 <- log2(tumor_Data[,"SYDE1"]+1)
CLIP3 <- log2(tumor_Data[,"CLIP3"]+1)
MRVI1 <- log2(tumor_Data[,"MRVI1"]+1)
PKN2 <- log2(tumor_Data[,"PKN2"]+1)
AHSA2 <- log2(tumor_Data[,"AHSA2"]+1)
AKAP11 <- log2(tumor_Data[,"AKAP11"]+1)
TIMP2 <- log2(tumor_Data[,"TIMP2"]+1)
CDK1 <- log2(tumor_Data[,"CDK1"]+1)
ABCE1 <- log2(tumor_Data[,"ABCE1"]+1)
PKD1 <- log2(tumor_Data[,"PKD1"]+1)
SGMS2 <- log2(tumor_Data[,"SGMS2"]+1)
MGP <- log2(tumor_Data[,"MGP"]+1)
HSPB8 <- log2(tumor_Data[,"HSPB8"]+1)
BOC <- log2(tumor_Data[,"BOC"]+1)

#DEG&positive
TMTC3 <- log2(tumor_Data[,"TMTC3"]+1)
FXYD6 <- log2(tumor_Data[,"FXYD6"]+1)
PDZD4 <- log2(tumor_Data[,"PDZD4"]+1)
SLC35A3 <- log2(tumor_Data[,"SLC35A3"]+1)
TMED7 <- log2(tumor_Data[,"TMED7"]+1)
SCAF1 <- log2(tumor_Data[,"SCAF1"]+1)
TUB <- log2(tumor_Data[,"TUB"]+1)
MYH11 <- log2(tumor_Data[,"MYH11"]+1)
C14orf132 <- log2(tumor_Data[,"C14orf132"]+1)
SPARCL1 <- log2(tumor_Data[,"SPARCL1"]+1)
TRO <- log2(tumor_Data[,"TRO"]+1)

#DEG&negative
C12orf48 <- log2(tumor_Data[,"C12orf48"]+1)
C14orf129 <- log2(tumor_Data[,"C14orf129"]+1)
C18orf32 <- log2(tumor_Data[,"C18orf32"]+1)
PDLIM7 <- log2(tumor_Data[,"PDLIM7"]+1)
COPS4 <- log2(tumor_Data[,"COPS4"]+1)
ADAMTSL3 <- log2(tumor_Data[,"ADAMTSL3"]+1)
FHL1 <- log2(tumor_Data[,"FHL1"]+1)
GPRASP1 <- log2(tumor_Data[,"GPRASP1"]+1)
HMCN1 <- log2(tumor_Data[,"HMCN1"]+1)
GBP4 <- log2(tumor_Data[,"GBP4"]+1)
JAK2 <- log2(tumor_Data[,"JAK2"]+1)
MXRA8 <- log2(tumor_Data[,"MXRA8"]+1)
SETD1A <- log2(tumor_Data[,"SETD1A"]+1)
RAB27B <- log2(tumor_Data[,"RAB27B"]+1)
TNRC6A <- log2(tumor_Data[,"TNRC6A"]+1)
NUMA1 <- log2(tumor_Data[,"NUMA1"]+1)
MRPL50 <- log2(tumor_Data[,"MRPL50"]+1)
ZNF24 <- log2(tumor_Data[,"ZNF24"]+1)
LONRF2 <- log2(tumor_Data[,"LONRF2"]+1)
ZNF767 <- log2(tumor_Data[,"ZNF767"]+1)
ARFIP1 <- log2(tumor_Data[,"ARFIP1"]+1)
USP33 <- log2(tumor_Data[,"USP33"]+1)
C5orf44 <- log2(tumor_Data[,"C5orf44"]+1)
ZNF720 <- log2(tumor_Data[,"ZNF720"]+1)
UBA3 <- log2(tumor_Data[,"UBA3"]+1)
LDB2 <- log2(tumor_Data[,"LDB2"]+1)
CDK10 <- log2(tumor_Data[,"CDK10"]+1)



temp_TF<-!is.na(Day.to.Death)
Observ<-Day.to.Death
Observ[]<-"Alive"                                 
Observ[temp_TF]<-"Death"
Observ<-as.factor(Observ)

Day.to.Death_0 <- Day.to.Death
Day.to.Death_0[!temp_TF] <- Day.to.follow.up[!temp_TF]


##SLC22A17
med1 <- median(SLC22A17, na.rm = TRUE)
SLC22A171<-SLC22A17

SLC22A171[SLC22A17<med1] <- "Low"
SLC22A171[SLC22A17>=med1] <- "High"

summary(SLC22A171)

file_name_plots=paste("C:/Users/korea/desktop/COAD/step3.Survival/Survival_Curve_by_SLC22A17.eps", sep="")
postscript(file=file_name_plots, paper="letter",horizontal =TRUE)
par(mfrow=c(1,1))
par(mar=c(4,3.9, 7, 1.3) )

fit.byrisk <- survfit( Surv(Day.to.Death_0, Observ=="Death") ~ SLC22A171 )
pvalue <- 1-pchisq(survdiff( Surv(Day.to.Death_0, Observ=="Death") ~ SLC22A171 )$chisq,2)

summary(fit.byrisk)

write.csv(fit.byrisk, "C:/Users/korea/desktop/COAD/step3.Survival/Survival_SLC22A17.csv", row.names=F)

plot(fit.byrisk ,conf.int=F, col=c("red","black","blue","green"),
     lty=c(1,2,1,2), 
     lwd=c(2,2,2,2),
     main=paste("SLC22A17  p = ",signif(pvalue,2),sep=""), 
     xlab="Days", ylab="P(Survive)",cex.main=2.5 ,cex.axis=1.2 , cex.lab=1.2  )

tmp<-table(SLC22A171)
tmp_2<-table(SLC22A171[Observ=="Death"])

level_list <- levels(SLC22A171)
median_list <- as.numeric(level_list)
for( level_index in 1:length(level_list) )
{
  median_list[level_index] <-round( median(Day.to.Death_0[level_list[level_index]==SLC22A171],na.rm=TRUE) ,1)
}

legend("topright",paste(names(tmp),": N=",tmp," (",median_list,")", sep=""), col=c("red","black","blue","green"), lty=c(1,2,1,2),lwd=c(2,2,2,2) )   

dev.off()


