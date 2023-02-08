library(GEOquery)
eSet <- getGEO("GSE120103",destdir = '.',getGPL = F) 
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
length(eSet)
b<- eSet[[1]]
raw_exprSet <- exprs(b)


#查看样本分组信息  平台类型##########################################################
pd =pData(b)#查看样本分组信息
#write.csv(pd,"GSE120103_cinical_data.csv")
###提取样本信息
library(stringr)
group_list =str_split(as.character(pd$source_name_ch1),'Group ',simplify =T)[,2]
table(group_list)
save(raw_exprSet,group_list,file='GSE120103_raw_expeSet.Rdata')
load("GSE120103_raw_expeSet.Rdata")

dim(raw_exprSet)
eSet[[1]]@annotation##---显示平台信息


#########ID转换#####################################################################
options('download.file.method.GEOquery'='auto')
options('GEOquery.inmemory.gpl'=FALSE)
library(GEOquery)
library(Biobase) 

#download GPL file, put it in the current directory, and load it:
gpl<-getGEO("GPL6480",destdir =PrimaryDirectory )

colnames(Table(gpl))
probe2symbol<-Table(gpl)[,c(1,7)]
names(probe2symbol)<-c("probe_id","symbol")

#一次性删除NA和空白
probe2symbol [probe2symbol ==""] <-NA
probe2symbol <-probe2symbol [complete.cases(probe2symbol),]
ids<-probe2symbol
plot(table(sort(table(ids$symbol))))#

###########数据的载入
load('GSE120103_raw_expeSet.Rdata')
#判断是否需要对数据取对数
ex <- raw_exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
exprSet=raw_exprSet  #

table(rownames(exprSet)%in%ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet)%in%ids$probe_id,] 
dim(exprSet) 

ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
dim(ids)
exprSet[1:5,1:5]

jimmy<-function(exprSet,ids){
  tmp=by(exprSet,
         ids$symbol,
         function(x)rownames(x)[which.max(rowMeans(x))])
  probes=as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet)%in%probes,]
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}
new_exprSet<-jimmy(exprSet,ids)
new_exprSet
print(dim(new_exprSet))
save(new_exprSet,group_list,file='GSE120103_new_expeSet.Rdata')


######校正，标准化##################################################################
rm(list=ls())
load(file='GSE120103_new_expeSet.Rdata')
exprSet=new_exprSet
exprSet[rownames(exprSet)=="GAPDH",]
boxplot(exprSet)

if(T){
  library(limma)
  normalized_expr = normalizeBetweenArrays(exprSet) # 分位数标准化 method默认为quantile
  p1 <- boxplot(normalized_expr,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F)
  title(main = list('Normalization',cex = 2 ,font = 2),
        xlab = list('Sample list',cex = 1.5,font = 2),
        ylab = '',line = 0.7)
  
  exprSet<-normalized_expr
  
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  
}
save(exprSet,group_list,file='GSE120103_for_dif_expeSet.Rdata')



######绘图##################################################################
rm(list=ls())
load(file='GSE120103_for_dif_expeSet.Rdata')
##画图展示单基因在不同组别的表达
###轻症、重症Mettl3
data_p = as.data.frame(t(exprSet))##行列转置    #注意这种转置方式：与使用melt函数之间区别
library(stringr)
data_p$group=str_split(as.character(rownames(data_p)),'_',simplify =T)[,1]

#查看lengend顺序
data_p$group <- factor(data_p$group,levels=c("2A)","2B)"))#调整坐标轴顺序


library(dplyr)
iso <- filter(data_p,group %in% c("2A)","2B)"))##提取某列中 的特定行

#Figure1C
head(ROC)
ROC<-data_p

roc1<-roc(ROC$group,ROC$METTL3);roc1

plot.roc(ROC$group,ROC$METTL3,
         main="Confidence intervals", percent=TRUE,
         ci=TRUE, # compute AUC (of AUC by default)
         print.auc=TRUE) # print the AUC (will contain the CI)


#Figure 1A
##color边框，fill 填充)
library(ggplot2)
#install.packages("ggsignif")
library(ggsignif)

compaired <- list(c("2A)","2B)"))

p <-ggplot(iso, aes(group,METTL3,fill=group)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  xlab("") +
  ylab(paste("Expression of ","METTL3"))+
  ggtitle("GSE120103 expression of METTL3") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    strip.background = element_blank(),
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5))+
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF")) +
  stat_compare_means(method = "t.test")
print(p)

#Figure S1A
rm(list=ls())
load(file='GSE120103_for_dif_expeSet.Rdata')
#-----多组箱式图
###筛选基因，通过merge 而不是which 如果是多条件的化
m6a_genelist<- read.csv("m6a_genelist_no METTL3.csv",header=TRUE,
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
head(m6a_genelist)
# 准备画图所需数据exp_L
library(reshape2)
head(exprSet)
exp_L = melt(exprSet)
head(exp_L)
colnames(exp_L)=c('symbol','sample','value')
head(exp_L)
library(dplyr)
library(stringr)
exp_L$group=str_split(as.character(exp_L$sample),'_',simplify =T)[,1]
iso <- filter(exp_L,group %in% c("2A)","2B)"))##提取某列中 的特定行
###提取数据
exp_met3<-merge(iso,m6a_genelist,by="symbol")
#查看lengend顺序
exp_met3$group <- factor(exp_met3$group,levels=c("2A)","2B)"))#调整坐标轴顺序
## 调整X轴分子的顺序   待整理
exp_met3$symbol <- factor(exp_met3$symbol,levels=c("METTL14","WTAP","RBM15","RBM15B",
                                                   "KIAA1429","CBLL1","ZC3H13","ALKBH5","FTO",
                                                   "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3",
                                                   "IGF2BP1","IGF2BP2","IGF2BP3","HNRNPA2B1",
                                                   "HNRNPC","FMR1","LRPPRC","ELAVL1"))
library(ggplot2)
require(cowplot)
require(tidyverse)
require(ggplot2)
require(ggsci)
require(ggpubr)
p<-ggplot(exp_met3, aes(x=symbol, y=value, fill = group)) +
  geom_boxplot(outlier.shape=21,lwd = 0.1,##置异常值的点的类型为21，默认的点类型是不能填充的
               outlier.size=1.5) + 
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF")) + 
  theme_bw() + 
  theme(#text=element_text(family="Times"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  stat_compare_means(aes(group = group), label = "p.signif")+
  ggtitle("GSE120103 expression of m6a genes")+
  guides(fill=FALSE)
print(p)
dev.off()#



