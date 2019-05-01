####highsort high
rm(list=ls())
#library(readxl)
library(network)
library(sqldf)

name_index <- c( "supp 1" )
#normal_lambda <- c( 0.31431 )


high_list <- as.data.frame(read.csv("E:/20180907/clinical/network/exist/high_low_list.csv"))
#info_un <- as.data.frame(read.csv("E:/20180907/network/exist/unduplicate_high_highsort/C5orf41.csv"))
info_du <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/duplicate_high_highsort/MSRB3.csv"))

file_name=paste("F:/20180907/clinical/network/exist/exist_in_high_lambda.txt")
result_table<-as.matrix(read.table(file_name, header=FALSE, sep=",")) 
file_name=paste("F:/20180907/clinical/network/exist/high_grade_240_neighbor.csv")
hub_gene <- as.data.frame(read.table(file_name, header=TRUE, sep=","))

hub_gene$degree <- as.numeric(hub_gene$degree)
#hub_gene <- hub_gene[c(order(-hub_gene$degree)),]
rownames(hub_gene) <- NULL
View(hub_gene[1:15,1])

gene_top <- hub_gene[10,1] #top 10
#gene_top
for(i in gene_top){
  #i
  result_table_c <- result_table[result_table[,1] %in% i | result_table[,2] %in% i,]  
  #result_table_c
  b <- as.vector(unique(rbind(as.matrix(unique(result_table_c[,1])), as.matrix(unique(result_table_c[,2])))))
  c <- length(b) 
  true_c<-matrix(rep(0, c*c), c, c)
  
  colnames(true_c) <- b
  rownames(true_c) <- b
  
  node_info <- as.data.frame(b)
  colnames(node_info) <- c("name")
  node_info <- sqldf('SELECT a.name, b.type
                     FROM node_info a
                     LEFT JOIN high_list b USING(name)')

  #edge_info_un <- as.data.frame(b)
  #colnames(edge_info_un) <- c("name")
  #edge_info_un <- sqldf('SELECT c.name, d.type 
                        #FROM edge_info_un c 
                        #LEFT JOIN info_un d USING(name)')
  
  edge_info_du <- as.data.frame(b)
  colnames(edge_info_du) <- c("name")
  edge_info_du <- sqldf('SELECT e.name, f.type 
                        FROM edge_info_du e 
                        LEFT JOIN info_du f USING(name)')
  
  for(j in 1:dim(result_table_c)[1])
  {
    true_c[result_table_c[j,1], result_table_c[j,2]] <- result_table_c[j,3]
  }
  
  true_X<-true_c
  true_X[true_X!=0]<-1
  
  d<-as.numeric(result_table_c[,3])
  
  
  aa1 <- network(true_X)
  
  
  file_name=paste("F:/20180907/clinical/network/exist/supp1_high_highsort_",i,".png", sep="")
  png(filename=file_name, height=1000, width=1000, bg="white")
  plot(aa1,displaylabels=TRUE,
       boxed.labels=FALSE,
       arrowhead.cex=0,
       #label.cex= ifelse(is.na(edge_info_du$type),1.3,
                         #ifelse(edge_info_du$type=="b",4,1.3)), #character expansion factor for label text,
       #label.lwd = ifelse(edge_info_du$type == "b", 10 , 0),
       pad=0.25,
       label.pos=10,
       edge.lwd=d*3.5,
       vertex.cex=2.5,
       #vertex.cex= ifelse(is.na(edge_info_un$type) ,2.5,
                         #ifelse(edge_info_un$type == "a",2.5,2.5)), #expansion factor for vertices
       #edge.col = ifelse(is.na(edge_info_du$type),"black",
                         #ifelse(edge_info_du$type=="b","dark red","black")),
       label.col = ifelse(is.na(edge_info_du$type),"black",
                            ifelse(edge_info_du$type == "b","dark red","black")),
       #vertex.enclose = ifelse(is.na(edge_info_un$type),FALSE,
                               #ifelse(edge_info_un$type=="a",TRUE,FALSE)),
       vertex.col =  ifelse(is.na(node_info$type),"white",
                            ifelse(node_info$type == "up", "light pink", "light green")),
             
       mode ="fruchtermanreingold",  # "kamadakawai",  #"circle", #"fruchtermanreingold", 
       cex.main=3)
  
  dev.off()
}

######highsort low
rm(list=ls())
#library(readxl)
library(network)
library(sqldf)

name_index <- c( "supp 1" )
#normal_lambda <- c( 0.31431 )


high_list <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/high_low_list.csv"))
info_du <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/duplicate_high_highsort/MSRB3.csv"))

file_name=paste("F:/20180907/clinical/network/exist/exist_in_low_lambda.txt")
result_table<-as.matrix(read.table(file_name, header=FALSE, sep=",")) 
file_name=paste("F:/20180907/clinical/network/exist/low_grade_sorthigh_240_neighbor.csv")
hub_gene <- as.data.frame(read.table(file_name, header=TRUE, sep=","))

hub_gene$degree <- as.numeric(hub_gene$degree)
#hub_gene <- hub_gene[c(order(-hub_gene$degree)),]
rownames(hub_gene) <- NULL
View(hub_gene[1:15,1])
gene_top <- hub_gene[10,1] #top 10
#gene_top
for(i in gene_top){
  #i
  result_table_c <- result_table[result_table[,1] %in% i | result_table[,2] %in% i,]  
  #result_table_c
  b <- as.vector(unique(rbind(as.matrix(unique(result_table_c[,1])), as.matrix(unique(result_table_c[,2])))))
  c <- length(b) 
  true_c<-matrix(rep(0, c*c), c, c)
  
  colnames(true_c) <- b
  rownames(true_c) <- b
  
  node_info <- as.data.frame(b)
  colnames(node_info) <- c("name")
  node_info <- sqldf('SELECT a.name, b.type
                     FROM node_info a
                     LEFT JOIN high_list b USING(name)')
  
  #edge_info_un <- as.data.frame(b)
  #colnames(edge_info_un) <- c("name")
  #edge_info_un <- sqldf('SELECT c.name, d.type 
  #FROM edge_info_un c 
  #LEFT JOIN info_un d USING(name)')
  
  edge_info_du <- as.data.frame(b)
  colnames(edge_info_du) <- c("name")
  edge_info_du <- sqldf('SELECT e.name, f.type 
                        FROM edge_info_du e 
                        LEFT JOIN info_du f USING(name)')
  
  for(j in 1:dim(result_table_c)[1])
  {
    true_c[result_table_c[j,1], result_table_c[j,2]] <- result_table_c[j,3]
  }
  
  true_X<-true_c
  true_X[true_X!=0]<-1
  
  d<-as.numeric(result_table_c[,3])
  
  
  aa1 <- network(true_X)
  
  
  file_name=paste("F:/20180907/clinical/network/exist/supp1_nega_highsort_",i,".png", sep="")
  png(filename=file_name, height=1000, width=1000, bg="white")
  plot(aa1,displaylabels=TRUE,
       boxed.labels=FALSE,
       arrowhead.cex=0,
       #label.cex= ifelse(is.na(edge_info_du$type),1.3,
       #ifelse(edge_info_du$type=="b",4,1.3)), #character expansion factor for label text,
       #label.lwd = ifelse(edge_info_du$type == "b", 10 , 0),
       pad=0.25,
       label.pos=10,
       edge.lwd=d*3.5,
       vertex.cex=2.5,
       #vertex.cex= ifelse(is.na(edge_info_un$type) ,2.5,
       #ifelse(edge_info_un$type == "a",2.5,2.5)), #expansion factor for vertices
       #edge.col = ifelse(is.na(edge_info_du$type),"black",
       #ifelse(edge_info_du$type=="b","dark red","black")),
       label.col = ifelse(is.na(edge_info_du$type),"black",
                          ifelse(edge_info_du$type == "b","dark red","black")),
       #vertex.enclose = ifelse(is.na(edge_info_un$type),FALSE,
       #ifelse(edge_info_un$type=="a",TRUE,FALSE)),
       vertex.col =  ifelse(is.na(node_info$type),"white",
                            ifelse(node_info$type == "up", "light green", "light pink")),
       
       mode ="fruchtermanreingold",  # "kamadakawai",  #"circle", #"fruchtermanreingold", 
       cex.main=3)
  
  dev.off()
}

######negasort nega
rm(list=ls())
#library(readxl)
library(network)
library(sqldf)

name_index <- c( "supp 1" )
#normal_lambda <- c( 0.31431 )


high_list <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/high_low_list.csv"))
#info_un <- as.data.frame(read.csv("F:/20180907/network/exist/AFF4.csv"))
info_du <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/duplicate_nega_negasort/XPO1.csv"))

file_name=paste("F:/20180907/clinical/network/exist/exist_in_low_lambda.txt")
result_table<-as.matrix(read.table(file_name, header=FALSE, sep=",")) 
file_name=paste("F:/20180907/clinical/network/exist/low_grade_353_neighbor.csv")
hub_gene <- as.data.frame(read.table(file_name, header=TRUE, sep=","))

hub_gene$degree <- as.numeric(hub_gene$degree)
#hub_gene <- hub_gene[c(order(-hub_gene$degree)),]
rownames(hub_gene) <- NULL
View(hub_gene[1:15,1])
gene_top <- hub_gene[10,1] #top 10
#gene_top
for(i in gene_top){
  #i
  result_table_c <- result_table[result_table[,1] %in% i | result_table[,2] %in% i,]  
  #result_table_c
  b <- as.vector(unique(rbind(as.matrix(unique(result_table_c[,1])), as.matrix(unique(result_table_c[,2])))))
  c <- length(b) 
  true_c<-matrix(rep(0, c*c), c, c)
  
  colnames(true_c) <- b
  rownames(true_c) <- b
  
  node_info <- as.data.frame(b)
  colnames(node_info) <- c("name")
  node_info <- sqldf('SELECT a.name, b.type
                     FROM node_info a
                     LEFT JOIN high_list b USING(name)')
  
  #edge_info_un <- as.data.frame(b)
  #colnames(edge_info_un) <- c("name")
  #edge_info_un <- sqldf('SELECT c.name, d.type 
  #FROM edge_info_un c 
  #LEFT JOIN info_un d USING(name)')
  
  edge_info_du <- as.data.frame(b)
  colnames(edge_info_du) <- c("name")
  edge_info_du <- sqldf('SELECT e.name, f.type 
                        FROM edge_info_du e 
                        LEFT JOIN info_du f USING(name)')
  
  for(j in 1:dim(result_table_c)[1])
  {
    true_c[result_table_c[j,1], result_table_c[j,2]] <- result_table_c[j,3]
  }
  
  true_X<-true_c
  true_X[true_X!=0]<-1
  
  d<-as.numeric(result_table_c[,3])
  
  
  aa1 <- network(true_X)
  
  
  file_name=paste("F:/20180907/clinical/network/exist/supp1_nega_negasort_",i,".png", sep="")
  png(filename=file_name, height=1000, width=1000, bg="white")
  plot(aa1,displaylabels=TRUE,
       boxed.labels=FALSE,
       arrowhead.cex=0,
       #label.cex= ifelse(is.na(edge_info_du$type),1.3,
       #ifelse(edge_info_du$type=="b",4,1.3)), #character expansion factor for label text,
       #label.lwd = ifelse(edge_info_du$type == "b", 10 , 0),
       pad=0.25,
       label.pos=10,
       edge.lwd=d*3.5,
       vertex.cex=2.5,
       #vertex.cex= ifelse(is.na(edge_info_un$type) ,2.5,
       #ifelse(edge_info_un$type == "a",2.5,2.5)), #expansion factor for vertices
       #edge.col = ifelse(is.na(edge_info_du$type),"black",
       #ifelse(edge_info_du$type=="b","dark red","black")),
       label.col = ifelse(is.na(edge_info_du$type),"black",
                          ifelse(edge_info_du$type == "b","dark red","black")),
       #vertex.enclose = ifelse(is.na(edge_info_un$type),FALSE,
       #ifelse(edge_info_un$type=="a",TRUE,FALSE)),
       vertex.col =  ifelse(is.na(node_info$type),"white",
                            ifelse(node_info$type == "up", "light green", "light pink")),
       
       mode ="fruchtermanreingold",  # "kamadakawai",  #"circle", #"fruchtermanreingold", 
       cex.main=3)
  
  dev.off()
}
#########negasort high

rm(list=ls())
#library(readxl)
library(network)
library(sqldf)

name_index <- c( "supp 1" )
#normal_lambda <- c( 0.31431 )


high_list <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/high_low_list.csv"))
info_du <- as.data.frame(read.csv("F:/20180907/clinical/network/exist/duplicate_nega_negasort/XPO1.csv"))

file_name=paste("F:/20180907/clinical/network/exist/exist_in_high_lambda.txt")
result_table<-as.matrix(read.table(file_name, header=FALSE, sep=",")) 
file_name=paste("F:/20180907/clinical/network/exist/high_grade_sortlow_353_neighbor.csv")
hub_gene <- as.data.frame(read.table(file_name, header=TRUE, sep=","))

hub_gene$degree <- as.numeric(hub_gene$degree)
#hub_gene <- hub_gene[c(order(-hub_gene$degree)),]
rownames(hub_gene) <- NULL
View(hub_gene[1:15,1])
gene_top <- hub_gene[10,1] #top 10
#gene_top
for(i in gene_top){
  #i
  result_table_c <- result_table[result_table[,1] %in% i | result_table[,2] %in% i,]  
  #result_table_c
  b <- as.vector(unique(rbind(as.matrix(unique(result_table_c[,1])), as.matrix(unique(result_table_c[,2])))))
  c <- length(b) 
  true_c<-matrix(rep(0, c*c), c, c)
  
  colnames(true_c) <- b
  rownames(true_c) <- b
  
  node_info <- as.data.frame(b)
  colnames(node_info) <- c("name")
  node_info <- sqldf('SELECT a.name, b.type
                     FROM node_info a
                     LEFT JOIN high_list b USING(name)')
  
  #edge_info_un <- as.data.frame(b)
  #colnames(edge_info_un) <- c("name")
  #edge_info_un <- sqldf('SELECT c.name, d.type 
  #FROM edge_info_un c 
  #LEFT JOIN info_un d USING(name)')
  
  edge_info_du <- as.data.frame(b)
  colnames(edge_info_du) <- c("name")
  edge_info_du <- sqldf('SELECT e.name, f.type 
                        FROM edge_info_du e 
                        LEFT JOIN info_du f USING(name)')
  
  for(j in 1:dim(result_table_c)[1])
  {
    true_c[result_table_c[j,1], result_table_c[j,2]] <- result_table_c[j,3]
  }
  
  true_X<-true_c
  true_X[true_X!=0]<-1
  
  d<-as.numeric(result_table_c[,3])
  
  
  aa1 <- network(true_X)
  
  
  file_name=paste("F:/20180907/clinical/network/exist/supp1_high_negasort_",i,".png", sep="")
  png(filename=file_name, height=1000, width=1000, bg="white")
  plot(aa1,displaylabels=TRUE,
       boxed.labels=FALSE,
       arrowhead.cex=0,
       #label.cex= ifelse(is.na(edge_info_du$type),1.3,
       #ifelse(edge_info_du$type=="b",4,1.3)), #character expansion factor for label text,
       #label.lwd = ifelse(edge_info_du$type == "b", 10 , 0),
       pad=0.25,
       label.pos=10,
       edge.lwd=d*3.5,
       vertex.cex=2.5,
       #vertex.cex= ifelse(is.na(edge_info_un$type) ,2.5,
       #ifelse(edge_info_un$type == "a",2.5,2.5)), #expansion factor for vertices
       #edge.col = ifelse(is.na(edge_info_du$type),"black",
       #ifelse(edge_info_du$type=="b","dark red","black")),
       label.col = ifelse(is.na(edge_info_du$type),"black",
                          ifelse(edge_info_du$type == "b","dark red","black")),
       #vertex.enclose = ifelse(is.na(edge_info_un$type),FALSE,
       #ifelse(edge_info_un$type=="a",TRUE,FALSE)),
       vertex.col =  ifelse(is.na(node_info$type),"white",
                            ifelse(node_info$type == "up", "light pink", "light green")),
       
       mode ="fruchtermanreingold",  # "kamadakawai",  #"circle", #"fruchtermanreingold", 
       cex.main=3)
  
  dev.off()
}

