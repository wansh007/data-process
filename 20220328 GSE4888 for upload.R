rm(list=ls())
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

library(GEOquery)
eSet <- getGEO("GSE4888", 
               destdir = '.', 
               getGPL = F) 
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
length(eSet)
b<- eSet[[1]]
raw_exprSet <- exprs(b)#

###提取样本分组信息
pd =pData(b)
head(pd)
library(stringr)
group_list<-as.character(pd$characteristics_ch1)
group_list

save(raw_exprSet,group_list,file='GSE4888_raw_expeSet.Rdata')
load("GSE4888_raw_expeSet.Rdata")
##进行查看下载的完整性，以及进行注释
dim(raw_exprSet)###完整性
eSet[[1]]@annotation####---显示平台信息

####ID转换######################################
library(AnnoProbe)
ids<-AnnoProbe::idmap('GPL570',type='bioc')
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol)))) 

###########数据的载入
load('GSE4888_raw_expeSet.Rdata')

#########################################################################################
#判断是否需要对数据取对数
ex <- raw_exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
########################################################################################
exprSet=log2(raw_exprSet) ##需要取对数，使用该代码 
table(rownames(exprSet)%in%ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet)%in%ids$probe_id,] #删除样本中，没有eg2probe 探针id号的
dim(exprSet) 

ids=ids[match(rownames(exprSet),ids$probe_id),]####将 eg2probe里面的gene按照exprSet中的顺序排列
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
save(new_exprSet,group_list,file='GSE4888_new_expeSet.Rdata')

######校正，标准化
rm(list=ls())
load(file='GSE4888_new_expeSet.Rdata')
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
save(exprSet,group_list,file='GSE4888_for_dif_expeSet.Rdata')

#绘图
######相关性Figure7A Figure7B Figure7C 
rm(list=ls())
load(file='GSE4888_for_dif_expeSet.Rdata')
library(dplyr)
data_p = as.data.frame(t(exprSet))##
boxplot(exprSet)

recepivity_gene <- c("METTL3","PGR","ELF3","MYC")  

library(ggstatsplot)             
library(ggside)
library(ggplot2)
library(ggpubr)
for(i in recepivity_gene){
  print(
    ggplot(data_p, aes_string(x ="METTL3",y =i))+ 
      geom_point()+ geom_smooth(method = 'lm', se = T, color = 'red')+
      theme_bw()+stat_cor(data=data_p, method = "pearson",alternative = c("less"))+theme_bw() + 
      theme(panel.grid=element_blank())
  )}




