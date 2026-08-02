#读取 OTU 丰度表
setwd("C:/Users/26355/Desktop/xxx/alpha")
library(vegan)
library(picante)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
input='feature-table.tsv'
df <- read.delim(input,header = T, row.names = 1,check.names = F)
Shannon<-diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson<-diversity(df, index = "simpson", MARGIN = 2, base = exp(1))
Richness <- specnumber(df,MARGIN = 2) #spe.rich ==sobs
index<-as.data.frame(cbind(Shannon,Simpson,Richness))

tdf<-t(df)
tdf<-ceiling(as.data.frame(t(df)))
obs_chao_ace<-t(estimateR(tdf))

obs_chao_ace<-obs_chao_ace[rownames(index),]
index$Chao<-obs_chao_ace[,2]
index$Ace<-obs_chao_ace[,4]
index$Sobs<-obs_chao_ace[,1]

index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df == 1) / colSums(df)
write.table(cbind(sample=c(rownames(index)),index), 'diversity.index.txt',row.names = F,sep = '\t',quote = F)

