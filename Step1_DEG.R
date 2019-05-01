rm(list=ls())

#load data for high and low group in COAD data set
high_colon <- read.csv("C:/Users/korea/Desktop/COAD/clinical/high_grade_colon_mRNA_standardize.csv", row.names = 1)
low_colon <- read.csv("C:/Users/korea/Desktop/COAD/clinical/low_grade_colon_mRNA_standardize.csv", row.names = 1)


wilcoxon_test.p <- vector()
wilcoxon_test.test <- vector()

## wilcoxon test
for (i in 1:ncol(high_colon)){
  wilcoxon_test.p[i] <- wilcox.test(high_colon[,i], low_colon[,i])$p.value
}


##DEg result
DEG <- colnames(high_colon)[wilcoxon_test.p<0.005]
DEG_result <- wilcoxon_test.p[wilcoxon_test.p<0.005]
total_DEG_result_0.005 <- cbind(DEG, DEG_result)
write.csv(total_DEG_result_0.005, "C:/Users/korea/Desktop/COAD/step1.DEG/DEG_result_0.05.csv", row.names=F)
