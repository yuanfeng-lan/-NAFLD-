####合并GEO数据集####
#网址：https://cloud.tencent.com/developer/article/1521695
setwd("TCGA.KIRC")
setwd("batch2")
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("sva")
library(sva)
library(tidyverse)


GSE48452 <- read.table("GSE48452.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE89632 <- read.table("GSE89632.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE164760 <- read.table("GSE164760.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

GSE48452group <- read.table("GSE48452-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE89632group <- read.table("GSE89632-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE164760group <- read.table("GSE164760-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

boxplot(GSE48452)
boxplot(GSE89632)
boxplot(GSE164760)


GSE48452 <- GSE48452[,rownames(GSE48452group)]
GSE89632 <- GSE89632[,rownames(GSE89632group)]
GSE164760 <- GSE164760[,rownames(GSE164760group)]


save(GSE36895,file="GSE36895.Rda")

save(GSE53757,file="GSE53757.Rda")

GSE48452$symbol <- rownames(GSE48452)
GSE89632$symbol <- rownames(GSE89632)
GSE164760$symbol <- rownames(GSE164760)




load("GSE36895.Rda")
load("GSE53757.Rda")

merge_eset=inner_join(GSE48452,GSE89632,by="symbol")

merge_eset=inner_join(merge_eset,GSE164760,by="symbol")

rownames(merge_eset) <- merge_eset$symbol
merge_eset <- merge_eset[,-74] #删除gene symble列
boxplot(merge_eset,outline=FALSE, notch=T,col=modType, las=2)
boxplot(merge_eset,outline=F)
dim(merge_eset)
exp <- as.matrix(merge_eset)
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(data)
class(data)
batchType <- c(rep(1,73),rep(2,63),rep(3,133))
modType <- c(rep("HC",14),rep("HO",27),rep("SS",14),rep("NASH",18),rep("HC",24),rep("SS",19),rep("NASH",20),rep("HC",6),rep("NASH",74),rep("HCC",53))
mod  <-  model.matrix(~as.factor(modType))
outTab <- data.frame(ComBat(data, batchType,mod, par.prior=TRUE))

modType = factor(modType,
                 levels = c("HC","HO",'SS','NASH','HCC'))

boxplot(outTab,outline=FALSE, notch=T,col=modType, las=2)

write.table(outTab, file = "normalize.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




####一个探针对应多个基因####
GSE164760$symbol <- data.frame(sapply(GSE164760$symbol,
                                    function(x)unlist(strsplit(x,'///'))[1]),
                             stringsAsFactors = F)[,1]
GSE164760 <- aggregate(GSE164760,)
GSE164760 <- aggregate(GSE164760[,1:170],by=list(GSE164760$symbol),mean,na.rm= TRUE)#排除最后的一列group
rownames(GSE164760) <- GSE164760$Group.1
GSE164760 <- GSE164760[,-1]
write.table(GSE164760, file = "GSE164760.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####PCA主成分分析####
install.packages("vegan")
library (vegan) #加载vegan包
library (ggplot2)#加载ggplot包

color=c( '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')
                    
otu <- read.delim("normalize.txt", sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
otu <- outTab
otu<-t(otu)
##无效
#otu <- exp0[,(1:149)]

otu<-as.matrix(otu)
otu<- ExpressionSet(assayData = otu)
#处理缺失值和异常值
otu <- filter.NA(otu, thres = 0.25)#排除超过25%的测量缺失的基因
otu <- fill.NA(otu, mode = 'mean')
otu <- filter.std(otu, min.std = 0)

otu <- as.data.frame(otu)
#对比表达谱数据 Hellinger 转化前后的差异
otu_hel <- decostand(otu, method = 'hellinger')
#使用 Hellinger 转化后的数据
pca_sp2 <- rda(otu_hel, scale = FALSE)

#特征值提取




pca_exp2 <- pca_sp2$CA$eig / sum(pca_sp2$CA$eig)
pc1 <- paste('PC1:', round(pca_exp2[1]*100, 2), '%')
pc2 <- paste('PC2:',round(pca_exp2[2]*100, 2), '%')
pca2=pca_sp2[["CA"]][["u"]][,c(1,2)]
map2<-read.table(file="group.txt",header=T,sep="\t",row.names=1)
group <- map2$group1
group <- map2$group2
merged2<-merge(pca2,map2,by="row.names",all.x=TRUE)
merged2$group1 = factor(merged2$group1,
                 levels = c("HC","HO",'SS','NASH','HCC'))
merged2$group2 = factor(merged2$group2,
                        levels = c("GSE48452","GSE164760",'GSE164760'))


p <- ggplot(data = merged2, aes(PC1, PC2)) +
  geom_point(aes(color = group1)) +  
  stat_ellipse(aes(fill =group1), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =color[1:length(unique(map2$group1))]) +
  scale_fill_manual(values = color[1:length(unique(map2$group1))]) +
  theme(panel.grid.major = element_line(color = 'gray', linewidth = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +  
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)
p

p <- ggplot(data = merged2, aes(PC1, PC2)) +
  geom_point(aes(color = group2)) +  
  stat_ellipse(aes(fill =group2), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =color[1:length(unique(map2$group2))]) +
  scale_fill_manual(values = color[1:length(unique(map2$group2))]) +
  theme(panel.grid.major = element_line(color = 'gray', linewidth = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +  
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)
p

####多组差异分析####

library(limma) 


data<-exp

group_list<-read.table(file="group.txt",header=T,sep="\t",row.names=1)
#new_group<- group_list[order(group_list[,1]),]  #对分组数据进行排序，按照数据框的第一列升序排序
group <- new_group#取第二列，这样处理得到的就是每一个位置所对应的样本属性信息

group <- group_list

suppressMessages(library(limma))


design <- model.matrix(~0+factor(map2$group1))
design <- model.matrix(~0+factor(group))

rownames(design) <- rownames(map2)

colnames(design)=c("HC","HCC",'HO','NASH','SS')


#design数据表要符合都是数字形式0、1，data看情况要转置

contrast.matrix<-makeContrasts("HCC-HC",levels=design)

##step1
fit <- lmFit(outTab,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 

write.table(nrDEG, file = "HCC-HC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(design, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####Mfuzz包时间趋势分析####
exp <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("genelightgreen.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

gene <- gene$gene
exp2 <- exp[,gene]


group <- group[rownames(exp),]

group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp2$group <- group$group1
dim(exp)
#数据处理
#[1]    16 12399
sample1<-aggregate(exp2[,1:62],by=list(exp2$group),mean,na.rm= TRUE)#排除最后的一列group

#更改顺序
sample1 <- sample1[c(1,3,5,4,2),]
#设置新行名
row.names(sample1)<-sample1[,1]
sample1<-data.frame(t(sample1[,-1]))

write.table(sample1, file = "samplelightgreen.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#安装R包
BiocManager::install("Mfuzz")
library("Mfuzz")
library(Mfuzz)
#构建对象
sample1<-as.matrix(sample1)
sample1<- ExpressionSet(assayData = sample1)

#处理缺失值和异常值
sample1 <- filter.NA(sample1, thres = 0.25)#排除超过25%的测量缺失的基因
sample1 <- fill.NA(sample1, mode = 'mean')
sample1 <- filter.std(sample1, min.std = 0)
#标准化
sample1 <- standardise(sample1)

#设置随机种子，设置需要展示的cluster的数量，然后聚类
set.seed(123)
cluster_num <- 6
sample1_cluster <- mfuzz(sample1, c = cluster_num, m = mestimate(sample1))

#作图
mfuzz.plot2(sample1, cl = sample1_cluster, mfrow = c(2, 3),
            time.labels = colnames(sample1),centre=TRUE,x11=F)
#导出基因
dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:10){
  potname<-names(sample1_cluster$cluster[unname(sample1_cluster$cluster)==i])
  write.csv(sample1_cluster[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}


####Mfuzz包时间趋势分析####
exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp0 <- exp[rownames(gene),]

group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp0 <- exp0 %>% t() %>% as.data.frame()

exp0$group <- group$group
dim(exp0)

#[1]    16 12399
sample1<-aggregate(exp0[,1:149],by=list(exp0$group),mean,na.rm= TRUE)#排除最后的一列group

#更改顺序
sample1 <- sample1[c(5,4,2,3,1),]
#设置新行名
row.names(sample1)<-sample1[,1]
sample1<-data.frame(t(sample1[,-1]))

#安装R包
BiocManager::install("Mfuzz")
library("Mfuzz")
library(Mfuzz)
#构建对象
sample1<-as.matrix(sample1)
sample1<- ExpressionSet(assayData = sample1)

#处理缺失值和异常值
sample1 <- filter.NA(sample1, thres = 0.25)#排除超过25%的测量缺失的基因
sample1 <- fill.NA(sample1, mode = 'mean')
sample1 <- filter.std(sample1, min.std = 0)
#标准化
sample1 <- standardise(sample1)

#设置随机种子，设置需要展示的cluster的数量，然后聚类
set.seed(123)
cluster_num <- 10
sample1_cluster <- mfuzz(sample1, c = cluster_num, m = mestimate(sample1))

#作图
mfuzz.plot2(sample1, cl = sample1_cluster, mfrow = c(2, 5),
            time.labels = colnames(sample1),centre=TRUE,x11=F)
#导出基因
dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:10){
  potname<-names(sample1_cluster$cluster[unname(sample1_cluster$cluster)==i])
  write.csv(sample1_cluster[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}
#最后，提取所有蛋白所属的聚类群，并和它们的原始表达值整合在一起
exp_cluster <- sample1_cluster$cluster
exp_cluster <- cbind(exp[names(exp_cluster), ], exp_cluster)
head(exp_cluster)
write.table(exp_cluster, file = "exp_cluster.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#如果您想提取数据分析过程中，标准化后的表达值（绘制曲线图用的那个值，而非原始蛋白表达值）
protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)
head(protein_standard_cluster)
#write.table(protein_standard_cluster, 'protein_standard_cluster.txt', sep = 't', col.names = NA, quote = FALSE)

####折线图的绘制####









####韦恩图绘制####


#从CRAN安装ggVennDiagram包；
install.packages("ggVennDiagram")
install.packages("ggsci")

#载入所需的R包；
library(ggplot2)
library(ggsci)
library(sf)
library(ggVennDiagram)
#'#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
#自定义颜色；
color1 <- alpha('#e6194B',0.9)
color2 <- alpha('#3cb44b',0.7)
color3 <- alpha('#ffe119',0.5)
color4 <- alpha('#4363d8',0.5)
color5 <- alpha('#f58231',0.5)
#绘制常见的4组数据venn图；
ggVennDiagram(x[1:4], label_alpha=0) +
  scale_fill_gradient(low="white",high =color4 ,guide="none")

#绘制5组数据的venn图；
#label_alpha = 0去除文字标签底色；
#category.names参数用于设定样本名称；
ggVennDiagram(x[1:5], label_alpha=0,label_size =3) +
  scale_color_brewer(palette = "Paired")+
  scale_fill_gradient(low="white",high = color1)

#edge_size设置线条粗细；
#label ="count"隐藏百分比，让原本就不“富裕”的空间更大一点；
#guide="none"隐藏图例；
ggVennDiagram(x[1:5], label_alpha=0,label_size =4,
              edge_size = 0.5,label ="count") +
  scale_color_brewer(palette = "Paired")+
  scale_fill_gradient(low="white",high = color1,guide="none")

#label = "none"用于隐藏图上的标签；
ggVennDiagram(x[1:5], label_alpha=0,label = "none",
              edge_size = 0.5) +
  scale_color_lancet()+
  scale_fill_gradient(low="gray100",high = "gray95",guide="none")

#还可以用交互的方式（plotly）查看每个子集中的基因；
ggVennDiagram(x[1:5], show_intersect = TRUE)

#下面再看下6组数据的绘制效果；
#个人感觉绘制6组时得到的图形区域有点胖，有点散，与文章的图还是有差异的；
ggVennDiagram(x[1:6], label_alpha=0,label_size =3,
              edge_size = 0.5,label ="count") +
  scale_color_lancet()+
  scale_fill_gradient(low="gray100",high = "gray95",guide="none")

ggVennDiagram(x, label_alpha=0,label_size =3,
              edge_size = 0.5,label ="count") +
  scale_color_lancet()+
  scale_fill_gradient(low="gray100",high = color2,guide="none")
#不添加过多的填充颜色，便于在Ai中进行后期调整；
ggVennDiagram(x, label_alpha=0,label = "none",
              edge_size = 0.5) +
  scale_color_lancet()+
  scale_fill_gradient(low="gray100",high = "gray95",guide="none")



####层次聚类####
datTraits <- traitData
datTraits <- as.numeric(unlist(datTraits))
#Re-cluster samples
sampleTree2= hclust(dist(dataExpr), method = "average")
#Convert traits to a color representation: white means low, red means high, greymeans missing entry
traitColors= numbers2colors(datTraits, signed = FALSE);
# Plotthe sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels= names(datTraits),
                    main ="Sample dendrogram and trait heatmap")
par(oma=c(8,8,9,8),mar=c(12,12,12,12))
plot(hh,xaxt='n')
par(cex=0.6)
axis(1)






####WGCNA2####
library(WGCNA)

#

library(reshape2)
library(stringr)

# 
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

## Allowing parallel execution with up to 47 working processes.

# 常规表达矩阵，log2转换后或
# Deseq2的varianceStabilizingTransformation转换的数据
# 如果有批次效应，需要事先移除，可使用removeBatchEffect
# 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，
# 需要quantile normalization

exprMat <- "WGCNA/LiverFemaleClean.txt"
exprMat <- 
  
  
exprMat <- read.table("normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp <- read.table("exp_p.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exprMat <- exprMat[rownames(exp),]

ids <- read.table("ids.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)




DEGs <- read.table("gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- DEGs$gene

group <- read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

sample <- subset(group,group!='SS') 
exprMat <- exprMat[,rownames(group)]
# 官方推荐 "signed" 或 "signed hybrid"
# 为与原文档一致，故未修改 
type = "unsigned"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
# 为与原文档一致，故未修改
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)
dataExpr <- exprMat
dim(dataExpr)

## [1] 3600  134
dataExpr <- dataExpr[gene,]
dataExpr <- subset(dataExpr,GSM2385720!='NA')


head(dataExpr)[,1:8]

##                 F2_2    F2_3     F2_14    F2_15    F2_19       F2_20
## MMT00000044 -0.01810  0.0642  6.44e-05 -0.05800  0.04830 -0.15197410
## MMT00000046 -0.07730 -0.0297  1.12e-01 -0.05890  0.04430 -0.09380000
## MMT00000051 -0.02260  0.0617 -1.29e-01  0.08710 -0.11500 -0.06502607
## MMT00000076 -0.00924 -0.1450  2.87e-02 -0.04390  0.00425 -0.23610000
## MMT00000080 -0.04870  0.0582 -4.83e-02 -0.03710  0.02510  0.08504274
## MMT00000102  0.17600 -0.1890 -6.50e-02 -0.00846 -0.00574 -0.01807182
##                F2_23    F2_24
## MMT00000044 -0.00129 -0.23600
## MMT00000046  0.09340  0.02690
## MMT00000051  0.00249 -0.10200
## MMT00000076 -0.06900  0.01440
## MMT00000080  0.04450  0.00167
## MMT00000102 -0.12500 -0.06820

## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))


dataExpr <- as.data.frame(t(dataExpr))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

## [1]  134 2697

head(dataExpr)[,1:8]

##       MMT00000051 MMT00000080 MMT00000102 MMT00000149 MMT00000159
## F2_2  -0.02260000 -0.04870000  0.17600000  0.07680000 -0.14800000
## F2_3   0.06170000  0.05820000 -0.18900000  0.18600000  0.17700000
## F2_14 -0.12900000 -0.04830000 -0.06500000  0.21400000 -0.13200000
## F2_15  0.08710000 -0.03710000 -0.00846000  0.12000000  0.10700000
## F2_19 -0.11500000  0.02510000 -0.00574000  0.02100000 -0.11900000
## F2_20 -0.06502607  0.08504274 -0.01807182  0.06222751 -0.05497686
##       MMT00000207 MMT00000212 MMT00000241
## F2_2   0.06870000  0.06090000 -0.01770000
## F2_3   0.10100000  0.05570000 -0.03690000
## F2_14  0.10900000  0.19100000 -0.15700000
## F2_15 -0.00858000 -0.12100000  0.06290000
## F2_19  0.10500000  0.05410000 -0.17300000
## F2_20 -0.02441415  0.06343181  0.06627665

#软阈值筛选
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

abline(h= 140, col = "red");#手动输入
#Determine cluster under the line
clust= cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)
#clust 1 contains the samples we want to keep.
keepSamples= (clust==1)
dataExpr= dataExpr[keepSamples, ]
nGenes= ncol(dataExpr)
nSamples= nrow(dataExpr)

write.table(dataExpr, file = "dataExpr.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

dataExpr <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#接下来我们引入临床数据，可视化展示临床的量化指标，同时综合以上的聚类图。此处用到一个重要的函数numbers2colors。
#此函数主要对数值化的参数进行高低的颜色标记，形成相应的热图。我们直接看下实例：
traitData= read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T);
traitData <- traitData[rownames(dataExpr),]

traitData <- traitData[,c(1,2,3,4,5,8,11,12,20,21,22)]

write.table(traitData, file = "traitData.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


dim(traitData)
names(traitData)
#remove columns that hold information we do not need.
allTraits= traitData[, -c(31, 16)];
allTraits= allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)

allTraits <- traitData[traitRows,]
# Forma data frame analogous to expression data that will hold the clinical traits.
femaleSamples= rownames(datExpr);
traitRows= match(femaleSamples, allTraits$Mice);
traitRows <- femaleSamples
datTraits= allTraits[traitRows,c(1:6)];
rownames(datTraits)= allTraits[traitRows, 1]


datTraits$gender <- ifelse(datTraits$gender=='female',1,0)


datTraits <- traitData
datTraits <- as.numeric(datTraits)
#Re-cluster samples
sampleTree2= hclust(dist(dataExpr), method = "average")
#Convert traits to a color representation: white means low, red means high, greymeans missing entry
traitColors= numbers2colors(datTraits, signed = FALSE);
# Plotthe sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels= names(datTraits),
                    main ="Sample dendrogram and trait heatmap")
par(oma=c(8,8,9,8),mar=c(12,12,12,12))
plot(hh,xaxt='n')
par(cex=0.6)
axis(1)



#计算软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)

## pickSoftThreshold: will use block size 2697.
##  pickSoftThreshold: calculating connectivity for given powers...
##    ..working on genes 1 through 2697 of 2697
##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.1370  0.825          0.412 587.000  5.95e+02  922.0
## 2      2   0.0416 -0.332          0.630 206.000  2.02e+02  443.0
## 3      3   0.2280 -0.747          0.920  91.500  8.43e+01  247.0
## 4      4   0.3910 -1.120          0.908  47.400  4.02e+01  154.0
## 5      5   0.7320 -1.230          0.958  27.400  2.14e+01  102.0
## 6      6   0.8810 -1.490          0.916  17.200  1.22e+01   83.7
## 7      7   0.8940 -1.640          0.869  11.600  7.29e+00   75.4
## 8      8   0.8620 -1.660          0.827   8.250  4.56e+00   69.2
## 9      9   0.8200 -1.600          0.810   6.160  2.97e+00   64.2
## 10    10   0.8390 -1.560          0.855   4.780  2.01e+00   60.1
## 11    12   0.8020 -1.410          0.866   3.160  9.61e-01   53.2
## 12    14   0.8470 -1.340          0.909   2.280  4.84e-01   47.7
## 13    16   0.8850 -1.250          0.932   1.750  2.64e-01   43.1
## 14    18   0.8830 -1.210          0.922   1.400  1.46e-01   39.1
## 15    20   0.9110 -1.180          0.926   1.150  8.35e-02   35.6
## 16    22   0.9160 -1.140          0.927   0.968  5.02e-02   32.6
## 17    24   0.9520 -1.120          0.961   0.828  2.89e-02   29.9
## 18    26   0.9520 -1.120          0.944   0.716  1.77e-02   27.5
## 19    28   0.9380 -1.120          0.922   0.626  1.08e-02   25.4
## 20    30   0.9620 -1.110          0.951   0.551  6.49e-03   23.5

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")



power = sft$powerEstimate
power

## [1] 6

# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}


##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 60,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

net <- blockwiseModules(
  dataExpr,
  power = 5,
  maxBlockSize = ncol(dataExpr),
  corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
  networkType = "unsigned",
  TOMType = "unsigned", 
  minModuleSize = 60,    ##越大模块越少
  mergeCutHeight = 0.2, ##越大模块越少
  numericLabels = TRUE, 
  saveTOMs= TRUE,
  saveTOMFileBase= "femaleMouseTOM",
  verbose = 3)

net <- recutBlockwiseTrees(dataExpr,
                           power = power,
                           maxBlockSize = ncol(dataExpr),
                           corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
                           networkType = "unsigned",
                           TOMType = "unsigned", 
                           minModuleSize = 30,    ##越大模块越少
                           mergeCutHeight = 0.25, ##越大模块越少
                           numericLabels = TRUE, 
                           saveTOMs= TRUE,
                           loadTOMs=TRUE,
                           saveTOMFileBase = paste0(exprMat, ".tom"),
                           verbose = 3)





##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will use 47 parallel threads.
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file WGCNA/LiverFemaleClean.txt.tom-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
##      ..removing 3 genes from module 1 because their KME is too low.
##      ..removing 5 genes from module 12 because their KME is too low.
##      ..removing 1 genes from module 14 because their KME is too low.
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules whose distance is less than 0.25
##        Calculating new MEs...

# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)

## 
##   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
## 135 472 356 333 307 303 177 158 102  94  69  66  63  62

#层级树展示各网络
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)


dataExpr <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exp <- dataExpr
exp <- as.data.frame(t(exp))
exp$gene <- rownames(exp)
exp$group <- moduleColors

colors_group <-exp[,c(269,270)] 

write.table(colors_group, file = "colors_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(colors_group, file = "colors_group2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


save(net, moduleColors, file = "step3_genes_modules.Rdata")
save(net, moduleColors, file = "step3_genes_modules2.Rdata")

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

#模块之间相关性图
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)


#TOM图
TOM=TOMsimilarityFromExpr(dataExpr,power=7)
dissTOM=1-TOM
## draw all genes 

geneTree = net$dendrograms[[1]]
plotTOM = dissTOM^7
diag(plotTOM)=NA
png("step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
TOMplot(plotTOM,geneTree,moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main="Network heapmap plot")

save(TOM, file = "TOM.Rdata")

# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main = "Network heatmap plot, all genes")

#有时基因太多，需要限制基因的数量
nSelect = 1000

nGenes=15026
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

TOMplot(plotDiss, selectTree, selectColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main="Network heapmap plot")


## 如果有表型数据，也可以跟ME数据放一起，一起出图
traitData <- 
  traitData <- read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
traitData <- traitData[rownames(dataExpr),]

traitData$`diagnosis:ch1` <- ifelse(traitData$`diagnosis:ch1`=='NASH',1,0)


MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)



# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")

save(TOM, file = "TOM.Rdata")

#导出用于cytoscape
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.5,
                               nodeNames = probes, nodeAttr = moduleColors)


#关联表型数据
trait <- read.table("phe2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trait <- trait[rownames(MEs),]
trait <- traitData
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
traitData <- trait
### 模块与表型数据关联
if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
win.graph(width=4.875, height=2.5,pointsize=8)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.3, 
               ySymbols = colnames(MEs_col), colorLabels = T, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.2, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## 从上图可以看到MEmagenta与Insulin_ug_l相关

## 模块内基因与表型数据关联

# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
# 值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要
# 。

### 计算模块与基因的相关性矩阵

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# 计算性状与基因的相关性矩阵

## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

####limma包geneSeT test计算模块内的基因和差异基因的相关性####

library(limma)

genetest <- read.csv("genesetTset.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
stat1 <- read.table("yellowset.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
k <- 1
stat <- rnorm(100)
sel <- 1:449; stat[sel] <- stat[sel]+1
sel <- genetest$Red
pval<-data.frame(t1=as.character(4:1)) 
for (i in c(1:4)) {
  stat <- stat1[,i]
  a[i] <- wilcoxGST(sel,stat)
  }


pval[k] <- a
k <- k+1

colnames(pval) <- c('Red','Black','Cyan','Blue','Brown','Green','Greenyellow','Magenta','Pink','Purple','Salmon','Tan','Turquoise','Yellow')
rownames(pval) <- c('HCC','NASH','SS','HO')
write.table(pval, file = "pval.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
pval2 <- -log10(pval)
write.table(pval2, file = "pval2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

stat <- stat1[,4]
stat <- as.matrix(stat)
wilcoxGST(sel,stat)


library(ggplot2)
pval2 <- t(pval2)
testData <- read.csv("bubbletestdata.csv")
ggplot(pval2, aes(x=number,y=sales,size=percent))  + 
  geom_point(colour = "steelblue") 

library(pheatmap)
pheatmap(pval2)

pval[1,] <- a
y <- matrix(rnorm(1000*6),1000,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

# First set of 20 genes are genuinely differentially expressed
index1 <- 1:20
y[index1,4:6] <- y[index1,4:6]+1

# Second set of 20 genes are not DE
index2 <- 21:40

camera(y, index1, design)
camera(y, index2, design)

camera(y, list(set1=index1,set2=index2), design, inter.gene.cor=NA)
camera(y, list(set1=index1,set2=index2), design, inter.gene.cor=0.01)


#Not run: 

download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p1.rdata", 
	      "human_c2_v5p1.rdata", mode = "wb")

load("human_c2_v5p1.rdata")
c2.indices <- ids2indices(Hs.c2, y$genes$GeneID)

camera(y, c2.indices, design)
 
#End(Not run)



# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "yellow"
  module = 'magenta'
    pheno = "diagnosis"
    modNames = substring(colnames(MEs_col), 3)
    # 获取关注的列
    module_column = match(module, modNames)
    pheno_column = match(pheno,colnames(traitData))
    # 获取模块内的基因
    moduleGenes = moduleColors == module
    
    sizeGrWindow(7, 7)
    par(mfrow = c(1,1))
    # 与性状高度相关的基因，也是与性状相关的模型的关键基因
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                       abs(geneTraitCor[moduleGenes, pheno_column]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", pheno),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    
    
    weight =as.data.frame(datTraits$diagnosis);
    names(weight) ="diagnosis"
    # 命名模块
    modNames =substring(names(MEs), 3)
    geneModuleMembership =as.data.frame(cor(dataExpr, MEs, use ="p"));
    MMPvalue =as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    names(geneModuleMembership) =paste("MM", modNames, sep="");
    names(MMPvalue) =paste("p.MM", modNames, sep="");
    geneTraitSignificance =as.data.frame(cor(dataExpr, weight, use ="p"));
    GSPvalue =as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) =paste("GS.",names(weight), sep="");
    names(GSPvalue) =paste("p.GS.",names(weight), sep="");
    
    
    
    
    
    names(dataExpr)
    # 返回brown模块中所有ID
    names(dataExpr)[moduleColors=="blue"]
    # 导入注释文件
    annot = read.csv(file = "GeneAnnotation.csv");
    dim(annot)
    names(annot)
    probes = names(datExpr)
    probes2annot = match(probes, annot$substanceBXH)
    # 统计没有注释到的基因
    sum(is.na(probes2annot))
    # 创建数据集，包含探测ID ，
    geneInfo0 = data.frame(substanceBXH = probes, 
                           geneSymbol = annot$gene_symbol[probes2annot], 
                           LocusLinkID = annot$LocusLinkID[probes2annot], 
                           moduleColor = moduleColors, 
                           geneTraitSignificance, GSPvalue)
    # 通过显著性对模块进行排序
    modOrder = order(-abs(cor(MEs, weight, use = "p")))
    # 添加模块成员
    for (mod in 1:ncol(geneModuleMembership))
    {
      oldNames = names(geneInfo0)
      geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, 
                                                             modOrder[mod]],
                             MMPvalue[, modOrder[mod]]);
      names(geneInfo0) = c(oldNames, paste("MM.", 
                                           modNames[modOrder[mod]], sep=""),
                           paste("p.MM.", modNames[modOrder[mod]], sep=""))
    }
    # 对基因进行排序
    geneOrder = order(geneInfo0$moduleColor, abs(geneInfo0$GS.weight))
    geneInfo = geneInfo0[geneOrder, ]
    # 导出
    write.csv(geneInfo, file = "geneInfo.csv")
    
    
    write.table(geneTraitCor, file = "geneTraitCor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    write.table(geneTraitP, file = "geneTraitP.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    
    
    
    
    
    
    
    
    
    
    
    
    ### 计算邻接矩阵
    adjacency = adjacency(dataExpr, power = power)
    
    ### 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM
    
    ### 层级聚类计算基因之间的距离树 
    geneTree = hclust(as.dist(dissTOM), method = "average")
    
    ### 模块合并
    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 30
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
    # Convert numeric lables into colors
    dynamicColors = labels2colors(dynamicMods)
    
    ### 通过计算模块的代表性模式和模块之间的定量相似性评估，合并表达图谱相似
    
    MEList = moduleEigengenes(datExpr, colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs)
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average")
    MEDissThres = 0.25
    
    # Call an automatic merging function
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    # Eigengenes of the new merged
    
    ## 分步法完结
    
    #提取指定模块的基因名
    
    # Recalculate topological overlap
    TOM = TOMsimilarityFromExpr(Expr, power = 6); 
    # Select module
    module = "grey";
    module = "blue";
    module = "greenyellow";
    probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
    inModule = (moduleColors==module);
    modProbes = probes[inModule];
    modTOM = TOM[inModule, inModule]
    dimnames(modTOM) = list(modProbes, modProbes)
    cyt = exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
      nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.05,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule]
    )
    

    ####功能富集分析结果####
    library(tidyverse)
    library(AnnotationDbi)
    library("BiocManager")
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(dplyr)
    library(ggplot2)
    
    
    
    hh <- read.table("KEGG.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    rownames(hh) <- 1:nrow(hh)
    hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
    ggplot(hh,aes(y=order,x=Count,fill=PValue))+
      geom_bar(stat = "identity",width=0.7)+####柱子宽度
      #coord_flip()+##颠倒横纵轴
      scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
      labs(
        x = "Gene numbers", 
        y = "Pathways")+
      theme(axis.title.x = element_text(face = "bold",size = 16),
            axis.title.y = element_text(face = "bold",size = 16),
            legend.title = element_text(face = "bold",size = 16))+
      theme_bw()
    ggplot(hh,aes(y=order,x=Count))+
      geom_point(aes(size=Count,color=-1*PValue))+# 修改点的大小
      scale_color_gradient(low="#66C3A5",high = "#FD8D62")+
      labs(color=expression(PValue,size="Count"), 
           x="Gene Number",y="Pathways")+
      theme_bw()
    
    hh <- read.table("GOBP.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    rownames(hh) <- 1:nrow(hh)
    hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
    ggplot(hh,aes(y=order,x=Count,fill=PValue))+
      geom_bar(stat = "identity",width=0.7)+####柱子宽度
      #coord_flip()+##颠倒横纵轴
      scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
      labs(
        x = "Gene numbers", 
        y = "GO term")+
      theme(axis.title.x = element_text(face = "bold",size = 16),
            axis.title.y = element_text(face = "bold",size = 16),
            legend.title = element_text(face = "bold",size = 16))+
      theme_bw()
    library(tidyverse)
    library("BiocManager")
    library(org.Hs.eg.db)
    library(clusterProfiler)
 
    gene <- read.table("gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    gene <- gene$gene
    genelist <- bitr(gene, fromType="SYMBOL",
                     toType="ENTREZID", OrgDb='org.Hs.eg.db')
    DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
    
    #GO分析
    ego <- enrichGO(gene = gene,
                    OrgDb = org.Hs.eg.db, 
                    ont = "all",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff =0.05, 
                    qvalueCutoff =0.05,
                    readable = TRUE)
    
    ego_res <- ego@result
    save(ego,ego_res,file = "GO_PDZK1_DEG.Rdata")
    
    #3. 可视化
    ##3.1 柱状图
    barplot(ego, showCategory = 20,color = "pvalue",col = '#66C3A5')
    ##3.2 气泡图
    dotplot(ego, showCategory = 20)
    ##3.3 分类展示
    barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scale='free')
    dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scale='free')
    
    
    ####PDZK1差异分析结果KEGG富集分析####
    setwd("PDZK1_KEGG")
    #install.packages("tidyverse")
    #install.packages("BiocManager")
    BiocManager::install('clusterProfiler')
    #BiocManager::install('org.Hs.eg.db')
    library(tidyverse)
    library("BiocManager")
    library(org.Hs.eg.db)
    library(clusterProfiler)
    DEG <- as.data.frame(res)%>% 
      arrange(padj) %>% 
      dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)
    
    DEG <- DEG %>% rownames_to_column("Gene")
    
    genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                     toType="ENTREZID", OrgDb='org.Hs.eg.db')
    DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
    
    #KEGG分析
    kk <- enrichKEGG(gene         = DEG$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = 0.1,
                     qvalueCutoff =0.1)
    kk_res <- kk@result
    save(kk,kk_res,file = "KEGG_PDZK1_DEG.Rdata")
    
    load("KEGG_PDZK1_DEG.Rdata")
    
    #柱状图
    barplot(kk, showCategory = 20,color = "pvalue")
    #气泡图
    dotplot(kk, showCategory = 20)
    
    dev.off()
    
    
    ####气泡图展示结果####
    BiocManager::install('clusterProfiler')
    
    
    
    library(tidyverse)
    library(AnnotationDbi)
    library("BiocManager")
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(dplyr)
    library(ggplot2)
    
    
    
    
    #准备好在DAVID上富集好的结果TXT文档
    ego <- read.table("GSE164760-GO-BP.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    DEGs <- read.table("gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    gene.df <- bitr(DEGs$gene,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
    
    
    #TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
    gene <- gene.df$ENTREZID
    
    
    ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                        OrgDb=org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "ALL",#富集的GO类型
                        pAdjustMethod = "BH",#这个不用管，一般都用的BH
                        minGSSize = 1,
                        pvalueCutoff = 0.05,#P值可以取0.05
                        qvalueCutoff = 0.05,
                        readable = TRUE)
    
    ego_CC <- enrichGO(gene = gene,
                       OrgDb=org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    
    ego_BP <- enrichGO(gene = gene,
                       OrgDb=org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    
    ego_MF <- enrichGO(gene = gene,
                       OrgDb=org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    
    #4、将结果保存到当前路径
    ego_ALL <- as.data.frame(ego_ALL)
    ego_result_BP <- as.data.frame(ego_BP)
    ego_result_CC <- as.data.frame(ego_CC)
    ego_result_MF <- as.data.frame(ego_MF)
    ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)#或者这样也能得到ego_ALL一样的结果
    write.csv(ego_ALL,file = "ego_ALL.csv",row.names = T)
    write.csv(ego_result_BP,file = "ego_result_BP.csv",row.names = T)
    write.csv(ego_result_CC,file = "ego_result_CC.csv",row.names = T)
    write.csv(ego_result_MF,file = "ego_result_MF.csv",row.names = T)
    write.csv(gene.df,file = "ego.csv",row.names = T)
    
    #5、但是有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，这时候跳过第4步，来到第5步
    display_number = c(30, 10, 10)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
    ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
    ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
    ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]
    
    ##将以上我们摘取的部分通路重新组合成数据框
    
    
    go_enrich_df <- data.frame(
      ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
      GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
      type=factor(c(rep("biological process", display_number[1]), 
                    rep("cellular component", display_number[2]),
                    rep("molecular function", display_number[3])), 
                  levels=c("biological process", "cellular component","molecular function" )))
    
    ##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
    for(i in 1:nrow(go_enrich_df)){
      description_splite=strsplit(go_enrich_df$Description[i],split = " ")
      description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
      go_enrich_df$Description[i]=description_collapse
      go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
    }
    
    ##开始绘制GO柱状图
    ###横着的柱状图
    go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
    COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色
    
    ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
      geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
      scale_fill_manual(values = COLS) + ###颜色
      coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
      xlab("GO term") + 
      ylab("Gene_Number") + 
      labs(title = "TCGA.LIHC-GO")+
      theme_bw()
    
    
    
    
    
    
    ###竖着的柱状图 
    go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
    COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
    ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
      geom_bar(stat="identity", width=0.8) + 
      scale_fill_manual(values = COLS) + 
      theme_bw() + 
      xlab("GO term") + 
      ylab("Num of Genes") + 
      labs(title = "The Most Enriched GO Terms")+ 
      theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置
    
    
    #1、KEGG富集
    kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
    
    #2、可视化
    ###柱状图
    hh <- as.data.frame(kk)#自己记得保存结果哈！
    
    
    
    hh <- read.table("KEGG.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    rownames(hh) <- 1:nrow(hh)
    hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
    ggplot(hh,aes(y=order,x=Count,fill=PValue))+
      geom_bar(stat = "identity",width=0.7)+####柱子宽度
      #coord_flip()+##颠倒横纵轴
      scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
      labs(title = "KEGG
       GSE167523",
           x = "Gene numbers", 
           y = "Pathways")+
      theme(axis.title.x = element_text(face = "bold",size = 16),
            axis.title.y = element_text(face = "bold",size = 16),
            legend.title = element_text(face = "bold",size = 16))+
      theme_bw()
    ggplot(hh,aes(y=order,x=Count))+
      geom_point(aes(size=Count,color=-1*PValue))+# 修改点的大小
      scale_color_gradient(low="#66C3A5",high = "#FD8D62")+
      labs(color=expression(PValue,size="Count"), 
           x="Gene Number",y="Pathways",title="GOBP
       GSE167523")+
      theme_bw()
    
    
####富集弦图####
BiocManager::install("GOplot")
    
#加载GOplot
library(GOplot)
    
    #加载测试数据
    data(EC)
    #创建circ对象，EC$david为富集分析结果
    #EC$genelist为差异表达分析结果
    GOBPdavid <- read.table("GOBPdavid.txt", sep = "\t",row.names = 1,check.names = F,header = T)
    GOBPdavid <- read.table("KEGGdavid.txt", sep = "\t",row.names = 1,check.names = F,header = T)
    genelist1 <- read.table("HCC-HC.txt", sep = "\t",row.names = 1,check.names = F,header = T)
   genelist1 <- genelist1[gene$gene,]
    genelist1$ID <- rownames(genelist1)
    genelist1 <- genelist1[,c(7,1,2,3,4,5,6)]
    
    
    gene <- read.table("gene.txt", sep = "\t",row.names = 1,check.names = F,header = T)
    
    circ <- circle_dat(GOBPdavid, genelist1)
    genes <- genelist1[,(1:2)]
    genes <- EC$genes
    genelist <- EC$genelist
    david <- EC$david
    circ <- circle_dat(david, genelist)
    process <- EC$process
    process <- GOBPdavid[,'Term']
    process <- process[1:20]
    #常见chord对象，EC$genes为需要展示的基因包含FC
    #EC$process为需要展示的GO条目
    chord <- chord_dat(circ, genes, process)
    
    #创建pdf文件，保存弦图
    pdf("chord_demo.pdf",height = 14,width = 13)
    #绘制弦图
    GOChord(chord,   #chord对象
            space = 0.02,  #右侧色块之间的间距
            gene.order = 'logFC',   #基因展示顺序根据logFC来
            gene.space = 0.25,  #基因名字和色块之间的距离
            gene.size = 5 #基因名字大小
           
    )
    dev.off()
    
    
####热图####
    library(pheatmap)
    
    gene <- read.table("gene4个.txt", sep = "\t",row.names = 1,check.names = F,header = T)
    gene <- read.table("genegreen.txt", sep = "\t",row.names = 1,check.names = F,header = T)
    gene <- gene$gene
    exp2 <-as.matrix(exp2)
    exp<- ExpressionSet(assayData = exp)
    exp_diff <- exp2
    exp_diff <- as.numeric(exp_diff[-102,])
    cg = rownames(DEG)[DEG$change !="NOT"]
    exp_diff <- exp[cg,]
    group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
    annotation_col=data.frame(group=group_list)
    rownames(annotation_col)=colnames(exp_diff)
    
    
    write.table(exp2, file = "PDZ_exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    exp_diff <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    exp_diff <- t(exp_diff)
    exp_diff2 <- exp_diff[gene,]
    sample <- colnames(exp_diff2)
    
    group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    group2 <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    group <- group[,1,drop=F]
    group <- group[colnames(exp_diff2),,drop=F]
    exp_diff2 <- exp_diff2[,rownames(group2)]
    write.table(group, file = "group1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    group$group1 <- as.factor(group$group1)
    group2$group1 = factor(group2$group1,
                     levels = c("HC","HO",'SS','NASH','HCC'))
    pheatmap(exp_diff2,
             annotation_col=group2,
             scale = "row",
             show_rownames = F,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =T,
             fontsize = 10,
             fontsize_row=3,
             fontsize_col=3)
    dev.off()
    
    
    
    GSE48452 <- read.table("GSE48452.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    GSE48452 <- GSE48452[gene,]
    GSE48452 <- GSE48452[,rownames(GSE48452group)]
    GSE48452group <- read.table("GSE48452-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    GSE48452group$group = factor(GSE48452group$group,
                                 levels = c("Control","Healthy_obese",'Steatosis','Nash'))
    
    pheatmap(GSE48452,
             annotation_col=GSE48452group,
             scale = "row",
             show_rownames = F,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =T,
             fontsize = 10,
             fontsize_row=3,
             fontsize_col=3)
    
    
    GSE89632 <- read.table("GSE89632.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    GSE89632 <- GSE164760[gene,]
   
    GSE89632group <- read.table("GSE89632-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    GSE89632 <- GSE89632[,rownames(GSE164760group)]
    GSE89632group$group = factor(GSE89632group$group,
                           levels = c("HC","SS",'NASH'))
    pheatmap(GSE89632,
             annotation_col=GSE89632group,
             scale = "row",
             show_rownames = F,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =T,
             fontsize = 10,
             fontsize_row=3,
             fontsize_col=3)
    
    GSE164760 <- read.table("GSE164760.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    GSE164760 <- GSE164760[gene,]
    
    GSE164760group <- read.table("GSE164760-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    GSE164760 <- GSE164760[,rownames(GSE164760group)]
    
    pheatmap(GSE164760,
             annotation_col=GSE164760group,
             scale = "row",
             show_rownames = F,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =T,
             fontsize = 10,
             fontsize_row=3,
             fontsize_col=3)
    
####多组样本t检验####
  exp3 <- exp[,c('PDZK1','MAGI3')]
  exp3$group <- group2$group1
  exp3$ID <- rownames(exp3)
  exp3 %>% sample_n_by(group, size = 1)
  
  exp3 <- exp3 %>%
    gather(key = "group", value = "exp",PDZK1,MAGI3) %>%
    convert_as_factor(ID, exp)
  set.seed(123)
  exp3 %>% sample_n_by(gene, exp, size = 1)
  
  
  
  stat.test <- exp3 %>%
    group_by(gene) %>%
    pairwise_t_test(
      exp ~ group, paired = TRUE, 
      p.adjust.method = "bonferroni"
    ) %>%
    select(-df, -statistic, -p) # Remove details
  stat.test
    
    
    
    
    
      library(tidyverse)
    library(rstatix)
    library(ggpubr)
    set.seed(123)
    data("anxiety", package = "datarium")
    anxiety %>% sample_n_by(group, size = 1)

    anxiety <- anxiety %>%
      gather(key = "time", value = "score", t1, t2, t3) %>%
      convert_as_factor(id, time)
    set.seed(123)
    anxiety %>% sample_n_by(group, time, size = 1)
    
    stat.test <- anxiety %>%
      group_by(group) %>%
      pairwise_t_test(
        score ~ time, paired = TRUE, 
        p.adjust.method = "bonferroni"
      ) %>%
      select(-df, -statistic, -p) # Remove details
    stat.test
    
    # 绘制图形
    bxp <- ggboxplot(
      anxiety, x = "group", y = "score",
      color = "time", palette = "jco"
    )
    # 添加显著性标记
    stat.test <- stat.test %>% add_xy_position(x = "group")
    bxp + stat_pvalue_manual(
      stat.test, label = "p.adj.signif", 
      step.increase = 0.08
    )
    
    
    bxp + stat_pvalue_manual(
      stat.test, label = "p.adj.signif", 
      step.increase = 0.08, hide.ns = TRUE, tip.length = 0
    )
    
    
    
    df <- data.frame(
      a = c(3.53, 4.59, 4.34, 2.66, 3.59, 3.13, 2.64, 2.56, 3.50, 3.25, 
            3.30, 4.04, 3.53, 3.56, 3.85, 4.07, 3.52, 3.93, 4.19, 2.96, 
            1.37, 3.93, 2.33, 2.98, 4.00, 3.55, 2.96, 4.30, 4.16, 2.59),
      b = c(2.42, 3.36, 4.32, 2.34, 2.68, 2.95, 1.56, 3.11, 1.81, 1.77, 
            1.98, 2.63, 2.86, 2.93, 2.17, 2.72, 2.65, 2.22, 2.90, 2.97, 
            2.36, 2.56, 2.52, 2.27, 2.98, 3.72, 2.80, 3.57, 4.02, 2.31),
      c = c(2.86, 2.28, 2.39, 2.28, 2.48, 2.28, 3.21, 2.23, 2.32, 2.68, 
            2.66, 2.32, 2.61, 3.64, 2.58, 3.65, 2.66, 3.68, 2.65, 3.02, 
            3.48, 2.42, 2.41, 2.66, 3.29, 2.70, 3.04, 2.81, 1.97, 1.68),
      d = c(0.89, 1.06, 1.08, 1.27, 1.63, 1.89, 1.19, 2.17, 2.28, 1.72, 
            1.98, 1.74, 2.16, 3.37, 2.97, 1.69, 0.94, 2.11, 2.81, 2.52, 
            1.31, 2.51, 1.88, 1.41, 3.19, 1.92, 2.47, 1.02, 2.10, 3.71)
    )
    library(tidyr)
    df2 <- df %>%
      pivot_longer(a:d, # 需要转换的变量
                   names_to = "group", # 转换后的新变量名称
                   values_to = "value") # 转换后的新变量值
    
    # 统计检验
    stat.test <- exp3 %>%
      anova_test(PDZK1 ~ group) %>% # 进行t检验
      add_significance() # 增加统计星号
    
    # 统计检验
    stat.test <- exp3 %>%
      t_test(MAGI3 ~ group) %>% # 进行t检验
      add_significance() # 增加统计星号
    stat.test # 输出结果
    
    ## 绘制带P值的箱型图
    bxp <- ggboxplot(exp3, x = "group", y = "PDZK1", 
                     color = "group",
                     palette = 'lancet',add = 'dotplot',
                     ggtheme = theme_bw(),legend='none')
 
    ## 设置统计检验的分组和坐标位置
    stat.test <- stat.test %>% 
      add_xy_position(x = "group")
    
    bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                             tip.length = 0.01)+ ## 手动添加p值
      stat_compare_means(method = 'anova') ## 添加整体anova的结果，这更科学。
   
    
     bxp <- ggboxplot(
      exp3, x = "group", y = "MAGI3",
      color = "group", palette = "jco"
    )
    # 添加显著性标记
    stat.test <- stat.test %>% add_xy_position(x = "group")
    bxp + stat_pvalue_manual(
      stat.test, label = "p.adj.signif", 
      step.increase = 0.08
    )
    
    
    
    ####GSVA以及ssGSEA分析####
    library(tidyverse)
    library(GSEABase)
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GSVA")
    BiocManager::install("GSEABase")
    
    library(GSVA)
    library(pheatmap)
    
    
    #读取表达矩阵设置PDZK1分组
    exp <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    exp <- read.table("exp_FPKM.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    exp <- read.table("exp_P.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    exp <- t(exp)
    exp <- as.data.frame(exp)
    exp_PDZK1 <- exp['PDZK1',]
    gene <- "PDZK1"#每次运行只改这个基因名
    med=median(as.numeric(exp[gene,]))
    exp <- as.data.frame(exp)
    conditions=data.frame(sample=colnames(exp),
                          group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
      column_to_rownames("sample")
    conditions <- conditions[order(conditions$group), ]
    
    
    
    write.table(conditions, file = "PDZK1-group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    #读取gene set
    geneSet <- read.table("genesets.v2023.1.Hs.gmt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    geneSet <- read.csv("genesets.v2023.1.Hs.gmt",header = F,sep = "\t") # 用EXCEL打开删除NA列第二列
    
    geneSet <- readLines("genesets.v2023.2.Hs.gmt")
    geneSet <- geneSet[,-c(1,3)]
    
    
    d <- '"E:\桌面文件转移\R\TCGAteech\肝脏疾病相关\GSE89632\GSVA\fatty"' #存放gmt文件的路径
    gmtfs <- list.files("fatty.genesets.v2023.1.Hs.gmt",pattern = 'symbols.gmt')  # 路径下所有结尾为symbols.gmt文件
    gmtfs
    kegg_list <- getGmt(file.path(d,gmtfs[1])) 
    go_list <- getGmt(file.path(d,gmtfs[2]))
    
    
    class(geneSet)
    
    
    geneSet <- geneSet %>%
      column_to_rownames("V2")%>%t()
    a <- geneSet
    a <- a[1:nrow(a),]
    set <- colnames(a)
    l <- list()
    #i <- "Activated CD8 T cell"
    for (i in set) {
      x <-  as.character(a[,i])
      x <- x[nchar(x)!=0]
      x <-  as.character(x)
      l[[i]] <-x
    }
    save(l,file = "./MYC.gene_set.Rdata")
    
    #开始ssGSEA分析
    
    class(exp)
    exp <- as.matrix(exp)
    
    ssgsea<- gsva(exp, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
    #GSVA
    gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)
    
    #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
    
    #热图可视化
    write.table(gs.exp, file = "gs.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    library(pheatmap)
    group3 <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    PDZK1_group <- read.table("PDZK1-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    PDZK1_group <- read.table("PDZK1_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    group3 <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    PDZK1_group <- as.matrix(PDZK1_group)
    
    sample <- rownames(PDZK1_group)
    sample <- rownames(group3)
    ssgsea <- ssgsea[,sample]
    gs.exp <- gs.exp[,sample]
    annotation_col <- group3
    
    annotation_col <- PDZK1_group
    
    rownames(annotation_col) <- colnames(res.ssgsea)
    pheatmap(ssgsea,scale = "row",
             show_colnames = F,
             # 不展示行名
             cluster_rows = F,
             # 不对行聚类
             cluster_cols = F,
             # 不对列聚类
             annotation_col = annotation_col,
             # 加注释
             cellwidth = 5,
             cellheight = 5,
             # 设置单元格的宽度和高度
             fontsize = 5
    )
    pheatmap(ssgsea,
             show_colnames = F,
             # 不展示行名
             cluster_rows = F,
             # 不对行聚类
             cluster_cols = F,
             # 不对列聚类
             annotation_col = annotation_col,
             # 加注释
             cellwidth = 5,
             cellheight = 5,
             # 设置单元格的宽度和高度
             fontsize = 5
    )
    pheatmap(gs.exp,scale = "row",
             show_colnames = F,
             # 不展示行名
             cluster_rows = T,
             # 不对行聚类
             cluster_cols = F,
             # 不对列聚类
             annotation_col = annotation_col,
             # 加注释
             cellwidth = 5,
             cellheight = 5,
             # 设置单元格的宽度和高度
             fontsize = 5
    )
    group3$group1 = factor(group3$group1,
                           levels = c("HC","HO",'SS','NASH','HCC'))
    pheatmap(gs.exp,
             annotation_col=group3,
             scale = "row",
             show_rownames = T,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =F,
             fontsize = 10,
             fontsize_row=3,
             fontsize_col=3)
    #### 进行limma差异处理 ####
    ##设定 实验组exp / 对照组ctr
    
    library(limma) 
    
    gl
    exp="primed"
    ctr="naive"
    
    design <- model.matrix(~0+factor(group3$group))
    colnames(design) <- levels(factor(group3))
    rownames(design) <- colnames(gsva_mat)
    contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"exp/ctrl"
                                     levels = design)
    
    
    
    
    design <- model.matrix(~0+factor(conditions$group))
    rownames(design) <- rownames(conditions)
    colnames(design)=c('low','high')
    
    
    group3 <- PDZK1_group
    group_list <- group3$group
    design <- model.matrix(~0+factor(group_list))
    row.names(design) <- rownames(group3)
    colnames(design)=c('high','low')
    contrast.matrix<-makeContrasts("high-low",levels=design)
    ##step1
    rownames(design) <- colnames(exp)
    fit <- lmFit(exp,design)
    
    fit <- lmFit(gs.exp,design)
    #data的列名要等于design的行名
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2) 
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput)
    write.table(nrDEG, file = "fatty.high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(nrDEG, file = "ribosome.high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(nrDEG, file = "high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    annotation_row <- nrDEG[,c('logFC','adj.P.Val')]
    
    fit1 <- lmFit(gsva_mat,design)                 #拟合模型
    fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
    efit <- eBayes(fit2)                         #修正
    
    summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
    tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
    degs <- na.omit(tempOutput) 
    write.csv(degs,"gsva_go_degs.results.csv")
    
    
    
    
    pheatmap(gs.exp,scale = "row",
             show_colnames = F,
             # 不展示行名
             cluster_rows = F,
             # 不对行聚类
             cluster_cols = F,
             # 不对列聚类
             annotation_col = annotation_col,
             annotation_row = annotation_row,
             # 加注释
             cellwidth = 5,
             cellheight = 5,
             # 设置单元格的宽度和高度
             fontsize = 5,
             main='ribosome.GSE48452'
    )
    #### 对GSVA的差异分析结果进行热图可视化 #### 
    ##设置筛选阈值
    padj_cutoff=0.05
    log2FC_cutoff=log2(2)
    
    keep <- rownames(nrDEG[nrDEG$adj.P.Val < padj_cutoff, ])
    length(keep)
    gs.exp <- gs.exp[,sample]
    dat <- gs.exp[keep[1:length(keep)],] #选取前50进行展示
    
    
    annotation_row <- annotation_row[rownames(dat),]
    
    pheatmap::pheatmap(dat, 
                       fontsize_row = 8,
                       height = 10,
                       width=18,
                       annotation_col = annotation_col,
                       annotation_row = annotation_row,
                       show_colnames = F,
                       cluster_cols = F,
                       show_rownames = T,
                       main='fatty.GSE167523')
    pheatmap::pheatmap(dat, 
                       fontsize_row = 8,
                       height = 10,
                       width=18,
                       annotation_col = annotation_col,
                       annotation_row = annotation_row,
                       show_colnames = F,
                       cluster_cols = F,
                       show_rownames = T,
                       main='ribosome.GSE162694')
    
    
    BiocManager::install("clusterProfiler", update = TRUE, ask = FALSE)
    #趋势探查
    library("Mfuzz")
    library(Mfuzz)
    group <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    gs.exp <- read.table("gs.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    gs.exp <- t(gs.exp)
    gs.exp <- as.data.frame(gs.exp)
    gs.exp <- gs.exp[rownames(group),]
    gs.exp$group <- group$group1
    dim(exp)
    #数据处理
    #[1]    16 12399
    sample1<-aggregate(gs.exp[,1:180],by=list(gs.exp$group),mean,na.rm= TRUE)#排除最后的一列group
    
    #更改顺序
    sample1 <- sample1[c(1,3,5,4,2),]
    #设置新行名
    row.names(sample1)<-sample1[,1]
    sample1<-data.frame(t(sample1[,-1]))
    
    #构建对象
    sample1<-as.matrix(sample1)
    sample1<- ExpressionSet(assayData = sample1)
    
    #处理缺失值和异常值
    sample1 <- filter.NA(sample1, thres = 0.25)#排除超过25%的测量缺失的基因
    sample1 <- fill.NA(sample1, mode = 'mean')
    sample1 <- filter.std(sample1, min.std = 0)
    #标准化
    sample1 <- standardise(sample1)
    
    #设置随机种子，设置需要展示的cluster的数量，然后聚类
    set.seed(123)
    cluster_num <- 10
    sample1_cluster <- mfuzz(sample1, c = cluster_num, m = mestimate(sample1))
    
    #作图
    mfuzz.plot2(sample1, cl = sample1_cluster, mfrow = c(2, 5),
                time.labels = colnames(sample1),centre=TRUE,x11=F)
    #导出基因
    dir.create(path="mfuzz",recursive = TRUE)
    for(i in 1:10){
      potname<-names(sample1_cluster$cluster[unname(sample1_cluster$cluster)==i])
      write.csv(sample1_cluster[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
    }
    
  ####LASSO####
    install.packages('glmnet')
    library(glmnet)
    graphics.off()  # clear all graphs
    rm(list = ls()) 
    
    # 示例数据准备
    N = 63 # 观测数
    p = 6  # 变量数
    
    exp <- read.table("expSS.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    exp <- t(exp)
    gene <- read.table("gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    gene <- gene$gene
    exp <- exp[gene,]
    phe <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    exp2 <- as.data.frame(exp[gene,])
    phe <- phe[colnames(exp),]
    phe2 <- phe[,1,drop=F]
  sample <- rownames(phe)
  write.table(exp, file = "expHC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  
    phe2 <- phe2[phe2$group1==c('HC','HO','SS'),]
exp <- exp[,sample]
    # X variable
    exp2 = as.matrix(exp2)
    x <- exp2
    # 计算标准化前的均值和标准差
    colMeans(exp)    # mean
    apply(exp,2,sd)  # standard deviation
    
    # 标准化
    x = scale(x,center = T,scale = T)
    
    # 计算标准化后的均值和标准差
    colMeans(x)    # mean
    apply(x,2,sd)  # standard deviation
    
    #——————————————-
    # Y variable
    #——————————————-
    beta = c( 0.15, -0.33,  0.25, -0.25, 0.05,rep(0, p/2-5), 
              -0.25,  0.12, -0.125, rep(0, p/2-3))
    
    y <- phe2
    # Y variable, standardized Y
    y = x%*%beta + rnorm(n, sd=0.5)
    y = scale(y)
    y <- t(phe2)
    y <- as.matrix(y)
    x <- t(x)
    
    # Model
    # 当lambda = 0.01
    lambda <- 0.01
    # lasso
    la.eq <- glmnet(x, y, lambda=lambda,
                    family='gaussian', 
                    intercept = F, alpha=1) 
    # 当alpha设置为0则为ridge回归，将alpha设置为0和1之间则为elastic net     
    # 系数结果 (lambda=0.01)
    la.eq$beta[,1]
    
    # Lasso筛选变量动态过程图
    set.seed(123)
    la.eq <- glmnet(x, y, family="gaussian", 
                    intercept = F, alpha=1) 
    # plot
    plot(la.eq,xvar = "lambda", label = T)
    # 也可以用下面的方法绘制
    #matplot(log(la.eq$lambda), t(la.eq$beta),
    #               type="l", main="Lasso", lwd=2)
    # Run cross-validation & select lambda
    #————————————————
    set.seed(123)
    mod_cv <- cv.glmnet(x=x, y=y, family="gaussian", # 默认nfolds = 10
                        intercept = F, alpha=1)
    
    plot(mod_cv) 
    
    # lambda.min : the λ at which the minimal MSE is achieved.
    
    # lambda.1se : the largest λ at which the MSE is within one standard error of the minimal MSE.
    print(paste(mod_cv$lambda.min,
                log(mod_cv$lambda.min)))
    print(paste(mod_cv$lambda.1se,
                log(mod_cv$lambda.1se)))
    
    # 这里我们以lambda.min为最优 λ
    best_lambda <- mod_cv$lambda.min
    best_lambda
    
    # 最终模型的系数估计
    #find coefficients of best model
    best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
    coef(best_model,s=best_lambda)
    best_model$beta[,1]
    save(best_model,file = 'best_model.Rdata')
    
    # 新观测
    #输入新的样本以及观测值
    new = matrix(rnorm(6), nrow=1, ncol=6) 
    #使用 lasso 回归模型预测
    predict(best_model, s = best_lambda, newx = new)
    
    y_predicted <- predict(best_model, s = best_lambda, newx = new)
    #find SST and SSE
    mean(y)
    y_predicted=0.628
    sst <- sum((y - mean(y))^2)
    sse <- sum((y_predicted - y)^2)
    
    #find R-Squared
    rsq <- 1 - sse/sst
    rsq
    
    #我们把这几个系数拿出来组成广义线性方程
    xy <- cbind(x ,y)
    xy <- as.data.frame(xy)
    mod<-glm(diagnosis~PDZK1+FADS2,family="gaussian",data =xy)
    
    summary(mod)
    
    
    
    
    set.seed(123)
    cv_fit <- cv.glmnet(x=x, y=y, family="gaussian", # 默认nfolds = 10
                        intercept = F, alpha=1)
    
    plot(cv_fit) 
    save(cv_fit,file='cv_fit.Rdata')
    model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.lse) )
    re=cbind(y ,lasso.prob)
    write.table(re, file = "re1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    #得到预测结果
    dat=as.data.frame(re[,1:2])
    colnames(dat)=c('event','prob')
    dat$event=as.factor(dat$event)#画图时需要factor
    library(ggpubr) 
    p <- ggboxplot(dat, x = "event", y = "prob",
                   color = "event", palette = "jco",
                   add = "jitter")
    #  Add p-value
    p + stat_compare_means()
    #得出预测结果
    set.seed(123)
    fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)#构建模型
    fit$beta
    #一倍SE内的更简洁的模型,是22个miRNA
    #fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    #head(fit$beta)# 这里是40个miRNA
    choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
    rownames(gene) <- gene$gene
    choose_gene1 <- gene[choose_gene,,drop=F]
    write.table(choose_gene1, file = "choose_gene.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
    re=cbind(y ,lasso.prob)
    
    write.table(re, file = "re.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    re1 <- read.table("re.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    #得到预测结果
    dat=as.data.frame(re[,1:2])
    colnames(dat)=c('event','prob')
    dat$event=as.factor(dat$event)#画图时需要factor
    library(ggpubr) 
    p <- ggboxplot(dat, x = "event", y = "prob",
                   color = "event", palette = "jco",
                   add = "jitter")
    #  Add p-value
    p + stat_compare_means()
    #得出预测结果
    
    fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)#构建模型
    fit$beta
    #一倍SE内的更简洁的模型,是22个miRNA
    #fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    #head(fit$beta)# 这里是40个miRNA
    choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
    
    
    ####LASSO2####
    library(glmnet)
    library(foreign)
    bc <- read.spss("E:/r/Breast cancer survival agec.sav",
                    use.value.labels=F, to.data.frame=T)
    bc <- na.omit(bc)
    
    y<-as.matrix(bc[,8])
    x<-as.matrix(bc[,c(2:7,9:11)])
    
    f1 = glmnet(x, y, family="binomial", nlambda=100, alpha=1) #这里alpha=1为LASSO回归，如果等于0就是岭回归
    #参数 family 规定了回归模型的类型：
    #family="gaussian" 适用于一维连续因变量（univariate）
    ##family="mgaussian" 适用于多维连续因变量（multivariate）
    #family="poisson" 适用于非负次数因变量（count）
    #family="binomial" 适用于二元离散因变量（binary）
    #family="multinomial" 适用于多元离散因变量（category）
    #我们这里结局指标是2分类变量，所以使用binomial
    print(f1)#把f1结果输出
    
    plot(f1, xvar="lambda", label=TRUE)
    #我们可以把数据集取一部分进行验证（这步不做也可以）
    predict(f1, newx=x[2:5,], type = "response")
    #然后通过glmnet自带函数进行交叉检验，并输出图形
    cvfit=cv.glmnet(x,y)
    
    plot(cvfit)
    
    #我们这个图中有两条虚线，一个是均方误差最小时的λ值，一个是距离均方误差最小时一个标准误的λ值，有点拗口没关系，我们只要知道它是多少就可以了
    
    cvfit$lambda.min#求出最小值
    
    cvfit$lambda.1se#求出最小值一个标准误的λ值
    
    l.coef2<-coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact = F)
    
    l.coef1<-coef(cvfit$glmnet.fit,s=cvfit$lambda.1se,exact = F)
    
    l.coef1
    
    l.coef2
    
    #我们把这几个系数拿出来组成广义线性方程
    
    mod<-glm(status~age+pathsize+lnpos+pr,family="binomial",data = bc)
    
    summary(mod)
    
    
    ####lasso3####
    
    install.packages('lars')
    rm(list=ls())
    options(stringsAsFactors = F)
    
    Rdata_dir='Rdata/'
    Figure_dir='figures/'
    # 加载上一步从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。
    load( file = 
            file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
    )
    dim(expr)
    dim(meta)
    # 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
    # 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息
    
    # 这里需要解析TCGA数据库的ID规律，来判断样本归类问题。
    group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')
    
    table(group_list)
    exprSet=na.omit(expr)
    
    table(group_list)
    exprSet=na.omit(expr)
    dim(exprSet)
    load(  file = 
             file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
    )
    dim(exprSet) ## remove the nomral
    head(phe)
    exprSet[1:4,1:4]
    head(colnames(exprSet))
    head(phe$ID)
    ## 必须保证生存资料和表达矩阵，两者一致
    all(substring(colnames(exprSet),1,12)==phe$ID)
    
    
    library(lars) 
    library(glmnet) 
    x=t(log2(exprSet+1))#归一化
    y=phe$event
    #用基因的表达情况预测生死
    model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)#拉手回归模型
    print(model_lasso)
    # 列%Dev代表了由模型解释的残差的比例，对于线性模型来说就是模型拟合的R^2(R-squred)。
    # 它在0和1之间，越接近1说明模型的表现越好，
    # 如果是0，说明模型的预测结果还不如直接把因变量的均值作为预测值来的有效。
    
    head(coef(model_lasso, s=c(model_lasso$lambda[29],0.009)))
    #找到最有意义的一个点纳入模型，就是说没有把所有基因都放入模型里面，只是找到了放入基因数使得模型最好的那个点
    #这个29是用眼睛看来的，我们可以用下面第三个代码划线来看
    
    plot(model_lasso, xvar = "norm", label = TRUE)
    
    plot(model_lasso,xvar = "lambda", label = T)
    
    cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
    #这里选取了1000个，精确些
    plot(cv_fit)
    #再用得到的最佳的位置去建模，lASSO可以防止过度拟合
    
    
    #lasso.prob是通过模型预测每个样本的stage值
    model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
    lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
    re=cbind(y ,lasso.prob)
    
    
    write.table(re, file = "lasso模型预测.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    #得到预测结果
    dat=as.data.frame(re[,1:2])
    colnames(dat)=c('event','prob')
    dat$event=as.factor(dat$event)#画图时需要factor
    library(ggpubr) 
    p <- ggboxplot(dat, x = "event", y = "prob",
                   color = "event", palette = "jco",
                   add = "jitter")
    #  Add p-value
    p + stat_compare_means()
    #得出预测结果
    
    fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)#构建模型
    fit$beta
    #一倍SE内的更简洁的模型,是22个miRNA
    #fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    #head(fit$beta)# 这里是40个miRNA
    choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
    #选出用于建模的那些基因
    length(choose_gene)
    myexpr=x[,choose_gene]
    mysurv=phe[,c("days","event")]
    mysurv=phe2
    mysurv$days[mysurv$days< 1] = 1 
    # 详细代码参见这个网站https://github.com/jeffwong/glmnet/blob/master/R/coxnet.R#
    fit <- glmnet( myexpr, Surv(mysurv$days,mysurv$event), 
                   family = "cox") 
    fit <- glmnet( myexpr, mysurv$diagnosis, 
                   family = "cox")
    #cox模型需要有时间和event数据
    #用包自带的函数画图
    plot(fit, xvar="lambda", label = TRUE)
    plot(fit, label = TRUE)
    ## 如果需要打印基因名，需要修改函数，这里不展开。
    
    library(pheatmap) 
    choose_matrix=expr[choose_gene,]
    choose_matrix[1:4,1:4]
    n=t(scale(t(log2(choose_matrix+1))))  #scale()函数去中心化和标准化
    #对每个探针的表达量进行去中心化和标准化
    n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
    n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
    n[1:4,1:4]
    
    ## http://www.bio-info-trainee.com/1980.html
    y <- t(y)
    y <- as.data.frame(y)
    y$diagnosis<- as.factor(y$diagnosis)
    annotation_col = as.data.frame(t(y))
    rownames(annotation_col)=colnames(expr)
    
    annotation_col <- ifelse(y$diagnosis=='0','HC','NASH')
    
    y$diagnosis <- annotation_col
    y <- order(y$diagnosis)
    
    
    pheatmap(exp2,show_colnames = F,annotation_col = y,
             filename = 'lasso_genes.heatmap.png')
    
    
    
    pheatmap(exp2,
             annotation_col = y,
             scale = "row",
             show_rownames = T,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =F,
             fontsize = 10,
             fontsize_row=10,
             fontsize_col=3)
    
    
    
    
    library(ggfortify)
    df=as.data.frame(t(choose_matrix))
    df$group=group_list
    png('lasso_genes.pca.png',res=120)
    autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
    dev.off()
    
    ## 也可以尝试其它主成分分析的R包，视频就不继续没完没了的讲解了。
    
    
    library("FactoMineR")
    library("factoextra")  
    ## 这里的PCA分析，被该R包包装成一个简单的函数，复杂的原理后面讲解。
    dat.pca <- PCA(t(choose_matrix), graph = FALSE) #'-'表示“非”
    fviz_pca_ind(dat.pca,repel =T,
                 geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
                 col.ind =  group_list, # color by groups 颜色组
                 # palette = c("#00AFBB", "#E7B800"),
                 addEllipses = TRUE, # Concentration ellipses 集中成椭圆
                 legend.title = "Groups"
    )
    
    
    
    ####AUC模型####
    install.packages('ROCR')
    library(ROCR)
    re1 <- as.data.frame(re1)
    re2 <- ifelse(re$group1=='2',1,0)
    
    xy <- x
    xy <- xy[rownames(re1),]
    xy <- as.data.frame(xy)
    
    re=cbind(re1 ,xy)
    re$NASH <- re2
    re <- re[re$diagnosis!='1',]
    xy <- xy[xy$diagnosis!='1',]
    re$PDZK1 <- xy$PDZK1
    re$FADS2 <- xy$FADS2
    
    
    
    pred_min1 <- prediction(re[,2], re[,1])
    auc_min1 = performance(pred_min1,"auc")@y.values[[1]]
    
    pred_min2 <- prediction(re[,4], re[,3])
    auc_min2 = performance(pred_min2,"auc")@y.values[[1]]
    
    #求得AUC值
    perf_min1 <- performance(pred_min1,"tpr","fpr")
    perf_min2 <- performance(pred_min2,"tpr","fpr")
    
    plot(perf_min1,colorize=FALSE, col="blue")
    plot(perf_min2,colorize=FALSE, col="red",add=T)
    #绘图
    lines(c(0,1),c(0,1),col = "gray", lty = 4 )
    # y=x
    text(0.8,0.2, labels = paste0("AUC = ",round(auc_min1,3)))
    text(0.8,0.3, labels = paste0("PDZK1 AUC = ",round(auc_min2,3)))
    
    # 加AUC值
    
    pred_min1 <- prediction(re[,9], re[,3])
    auc_min1 = performance(pred_min1,"auc")@y.values[[1]]
    
    perf_min1 <- performance(pred_min1,"tpr","fpr")
    
    plot(perf_min1,colorize=FALSE, col="blue")
    
    lines(c(0,1),c(0,1),col = "gray", lty = 4 )
    
    text(0.8,0.2, labels = paste0("AUC = ",round(auc_min1,3)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ####其他数据集验证####
    new <- read.table("GSE66676.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    group <- read.table("GSE66676group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    exp <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    
    new <- read.table("GSE126848new.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    group <- read.table("re1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    group0 <- group
    group <- group0[,5,drop=F]
    new$ZNF878 <- c('1')
    colnames(new) <- gene
    
    write.table(new, file = "GSE66676new.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(group, file = "group1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    group <- read.table("group1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    group <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    gene <- read.table("gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    gene <- gene[,1]
    
    new <- read.table("GSE66676new.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    new <- exp2
    group <- subset(group,'group'== c('NASH','healthy'))
    group$type <- ifelse(group$group=='NASH',2,0)
    
    group <- group[,2,drop=F]
    group <- as.matrix(group)
    new <- new[,rownames(group)]
    new <- new[gene,]
    new <- as.matrix(new)
    new <- t(new)
    new <- as.data.frame(new)
    new = matrix(rnorm(62), nrow=1, ncol=62) 
    #使用 lasso 回归模型预测
    predict(best_model, s = best_lambda, newx = new)
    
    y_predicted <- predict(best_model, s = best_lambda, newx = new)
    #find SST and SSE
    mean(group)
    y_predicted=1.105889
    sst <- sum((group - mean(group))^2)
    sse <- sum((y_predicted - y)^2)
    
    #find R-Squared
    rsq <- 1 - sse/sst
    rsq
    
    
    
    
    
    new <- as.matrix(new)
    
    model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
    predict
    lasso.prob <- predict(cv_fit, newx=new , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
    lasso.prob <- lasso.prob[rownames(group),]
    re <- cbind(group,lasso.prob)
    
    write.table(re, file = "GSE66676re.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    #得到预测结果
    
    dat <- dat[rownames(group),]
    dat=as.data.frame(re[,1:2])
    colnames(dat)=c('disease subtype','prob')
    dat$`disease subtype`=as.factor(dat$`disease subtype`)#画图时需要factor
    library(ggpubr) 
    p <- ggboxplot(dat, x = "disease subtype", y = "prob",
                   color = "disease subtype", palette = "jco",
                   add = "jitter")+ 
      geom_boxplot()
    #  Add p-value
    p + stat_compare_means()
    #得出预测结果
    set.seed(123)
    fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)#构建模型
    fit$beta
    #一倍SE内的更简洁的模型,是22个miRNA
    #fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
    #head(fit$beta)# 这里是40个miRNA
    choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
    #选出用于建模的那些基因
    length(choose_gene)
    myexpr=x[,choose_gene]
    mysurv=phe[,c("days","event")]
    mysurv=phe2
    mysurv$days[mysurv$days< 1] = 1 
    # 详细代码参见这个网站https://github.com/jeffwong/glmnet/blob/master/R/coxnet.R#
    fit <- glmnet( myexpr, Surv(mysurv$days,mysurv$event), 
                   family = "cox") 
    fit <- glmnet( myexpr, mysurv$diagnosis, 
                   family = "cox")
    #cox模型需要有时间和event数据
    #用包自带的函数画图
    plot(fit, xvar="lambda", label = TRUE)
    plot(fit, label = TRUE)
    ## 如果需要打印基因名，需要修改函数，这里不展开。
    
    library(pheatmap) 
    choose_matrix=expr[choose_gene,]
    choose_matrix[1:4,1:4]
    n=t(scale(t(log2(choose_matrix+1))))  #scale()函数去中心化和标准化
    #对每个探针的表达量进行去中心化和标准化
    n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
    n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
    n[1:4,1:4]
    
    ## http://www.bio-info-trainee.com/1980.html
    y <- t(y)
    y <- as.data.frame(y)
    y$diagnosis<- as.factor(y$diagnosis)
    annotation_col = as.data.frame(t(y))
    rownames(annotation_col)=colnames(expr)
    
    annotation_col <- ifelse(y$diagnosis=='0','HC','NASH')
    
    y$diagnosis <- annotation_col
    
    
    
    pheatmap(exp2,show_colnames = F,annotation_col = y,
             filename = 'lasso_genes.heatmap.png')
    
    
    new <- t(new)
    pheatmap(new,
             annotation_col = group,
             scale = "row",
             show_rownames = T,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =F,
             fontsize = 10,
             fontsize_row=10,
             fontsize_col=3)
    
    
    
    
    library(ggfortify)
    df=as.data.frame(t(choose_matrix))
    df$group=group_list
    png('lasso_genes.pca.png',res=120)
    autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
    dev.off()
    
    
    install.packages('ROCR')
    library(ROCR)
    pred_min <- prediction(re[,2], re[,1])
    auc_min = performance(pred_min,"auc")@y.values[[1]]
    #求得AUC值
    perf_min <- performance(pred_min,"tpr","fpr")
    plot(perf_min,colorize=FALSE, col="blue") 
    #绘图
    lines(c(0,1),c(0,1),col = "gray", lty = 4 )
    # y=x
    text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
    # 加AUC值
####TCGA.LIHC赋值####
    exp <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    gene <- read.table("gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    gene <- gene$gene
    exp2 <- exp[gene,]
    write.table(exp2, file = "exp2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    exp2 <- read.table("exp2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    exp2 <- t(exp2)
    exp2 <- as.matrix(exp2)
    
    
    lasso.prob <- predict(best_model, newx=exp2 , s=best_lambda)
    write.table(lasso.prob, file = "lassoprobe2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

    ####生存分析####
    
    setwd("survival")
    surv <- exp_sur
    surv <- read.table("survivalgroup3.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    surv$OS.time <- surv$OS.time/30
    surv$OS <- as.numeric(unlist(surv$OS))#OS和OS.time数据numeric化，可数化。
    surv$PDZK1 <- as.numeric(unlist(surv$PDZK1))#OS和OS.time数据numeric化，可数化。
    
    #median或者平均数mean
    median(surv$altered)
    surv$group <- ifelse(surv$altered > median(surv$altered),"altered High","altered Low")
    surv$group <- factor(surv$group, levels = c("altered Low","altered High")) 
    class(surv$group)
    table(surv$group)
    
    
    median(surv$PDZK1)
    surv$group <- ifelse(surv$PDZK1 > median(surv$PDZK1),"PDZK1 High","PDZK1 Low")
    surv$group <- factor(surv$group, levels = c("PDZK1 Low","PDZK1 High")) 
    class(surv$group)
    table(surv$group)
    
    
    
    
    install.packages("survival")
    library(survival)
    fitd <- survdiff(Surv(OS.time, OS) ~ risk,
                     data      = surv,
                     na.action = na.exclude)
    pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
    
    #2.2 拟合生存曲线
    fit <- survfit(Surv(OS.time, OS)~ risk, data = surv)
    summary(fit)
    
    #3. 绘制生存曲线
    #方法1
    ###3. 设置颜色，坐标
    plot(fit, conf.int = T,
         col = c("blue", "red"),
         lwd = 2,
         xlab = "Time(Months)",
         ylab = "Survival probablity(%)"
    )
    ###添加标签
    legend("topright",
           title = "Group",
           c("Low", "High"),
           lwd = 2, lty = 1,
           col = c("blue", "red"))
    ###添加P值
    p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
    text(25, 0.1, p.lab)
    #前两个参数表示粘贴的标签横纵坐标位置
    dev.off()
    
    library(ROCR)
    pred_min <- prediction(surv[,1], surv[,3])
    auc_min = performance(pred_min,"auc")@y.values[[1]]
    #求得AUC值
    perf_min <- performance(pred_min,"tpr","fpr")
    plot(perf_min,colorize=FALSE, col="blue") 
    #绘图
    lines(c(0,1),c(0,1),col = "gray", lty = 4 )
    # y=x
    text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
    # 加AUC值
    
    
    #得到预测结果
    
    
    re <- read.table("TCGA.fibro数据2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    
    dat=as.data.frame(re[,1:2])
    colnames(dat)=c('prob','event')
    dat$event=as.factor(dat$event)#画图时需要factor
    library(ggpubr) 
    p <- ggboxplot(dat, x = "event", y = "prob",
                   color = "event", palette = "jco",
                   add = "jitter")
    #  Add p-value
    p + stat_compare_means()
    #得出预测结果
    ####方法2做生存分析图####
    #install.packages("survminer")
    library(survminer)
    ggsurvplot(fit,
               data = surv,
               pval = p.lab,
               conf.int = TRUE, # 显示置信区间
               #risk.table = TRUE, # 显示风险表
               risk.table.col = "strata",
               palette = "jco", # 配色采用jco
               legend.labs = c("Low", "High"), # 图例
               size = 1,
               xlim = c(0,120), # x轴长度，一般为0-10年
               break.time.by = 20, # x轴步长为20个月
               legend.title = "",
               surv.median.line = "hv", # 限制垂直和水平的中位生存
               ylab = "Survival probability (%)", # 修改y轴标签
               xlab = "Time (Months)", # 修改x轴标签
               #ncensor.plot = TRUE, # 显示删失图块
               ncensor.plot.height = 0.25,
               risk.table.y.text = FALSE)
    
    dev.off()
    
    ####基因之间的相关性分析####
   
    ####多个性状的相关性分析及可视化
    
    install.packages("Hmisc")
    install.packages("PerformanceAnalytics")
    dd = as.data.frame(matrix(rnorm(1000),100,10))
    
    exp0 <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    gene <- read.table("choose_gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    exp1 <- exp0[rownames(group),gene$gene]
    
    head(dd)
    
    # 计算相关系数及显著性
    
    library(Hmisc)#加载包
    res2 <- rcorr(as.matrix(exp1))
    res2
    
    # 可视化
    library(PerformanceAnalytics)#加载包
    chart.Correlation(exp1, histogram=T, pch=19)
    chart
    library(corrplot)
    #绘制一个下三角的热图，这个包的使用在之前的博客写过，这里一笔带过
    exp1 <- t(exp1)
    cor_matr = cor(exp1)
    cor_matr
    corrplot(cor_matr, type="upper", order="hclust", tl.col="black", tl.srt=30,addCoef.col = "grey60")
    
    
    
    rownames(cor_matr) = paste(1:14, colnames(cor_matr),sep = ' ')
    colnames(cor_matr) = as.character(1:14) # 设置矩阵的列变量名为数字
    
    col3 <- colorRampPalette(c('DodgerBlue3','white', "OrangeRed")) # 返回一个函数col3, 用来设置相关矩阵图的colorbar的分段
    corrplot(cor_matr, method = 'circle', diag = F, type = 'full', outline = F,
             col = col3(20), cl.lim = c(-1,1),addgrid.col = NA,
             tl.pos = 'lb',tl.cex = 0.75, tl.col = 'black', tl.srt = 0, tl.offset = 0.5) # 绘制相关矩阵
    axis(1,at = 1:14, labels = NA, pos = 14.5, tck = -0.01) # 在图片上边添加坐标轴，设置其刻度为位置，标签，画的位置，刻度的朝向
    axis(4,at = 1:14, labels = NA, pos = 0.5, tck = -0.01) # 在图片左边添加坐标轴，设置参数同上
    
    
    
    
    
    #热图
    library(pheatmap)
    exp2 <-as.matrix(exp1)
    exp<- ExpressionSet(assayData = exp)
    exp_diff <- exp2
    exp_diff <- as.numeric(exp_diff[-102,])
    cg = rownames(DEG)[DEG$change !="NOT"]
    exp_diff <- exp[cg,]
    group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
    annotation_col=data.frame(group=group_list)
    rownames(annotation_col)=colnames(exp_diff)
    
    
    write.table(exp2, file = "PDZ_exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
    
    exp_diff <- read.table("PDZ_exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    re <- read.table("re1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
    exp1 <- exp1[rownames(re),]
   
    re$group = factor(re$group,
                            levels = c('SS','NASH','HCC'))
    
    exp1 <- t(exp1)
    pheatmap(exp1,
             annotation_col=re,
             scale = "row",
             show_rownames = T,
             show_colnames =F,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             cluster_cols =F,
             fontsize = 10,
             fontsize_row=6,
             fontsize_col=5)
    dev.off()
    exp_diff <- t(exp_diff)
    pheatmap(exp_diff)
    
    
    
    
    
    
    ####相关性点图####
    library(ggplot2)
    library(ggpubr)
    library(ggpmisc)
    theme_set(ggpubr::theme_pubr()+
                theme(legend.position = "top"))
    
    # Load data
    data("mtcars")
    #mtcars是示例数据
    df <- mtcars
    # Convert cyl as a grouping variable
    df$cyl <- as.factor(df$cyl)
    # Inspect the data
    head(df[, c("wt", "mpg", "cyl", "qsec")], 4)
    
    
    b <- ggplot(df, aes(x = wt, y = mpg))
    # Scatter plot with regression line
    b + geom_point()+
      geom_smooth(method = "lm", color = "black", fill = "lightgray") 
    
    # Add a loess smoothed fit curve
    b + geom_point()+
      geom_smooth(method = "loess", color = "black", fill = "lightgray")
    
    b + geom_point(shape = 17)+
      geom_smooth(method = "lm", color = "black", fill = "lightgray")
    
    # Add regression line and confidence interval
    # Add correlation coefficient: stat_cor()
    ggscatter(df, x = "wt", y = "mpg",
              add = "reg.line", conf.int = TRUE,    
              add.params = list(fill = "lightgray")
              
    )+
      stat_cor(method = "pearson", 
               label.x = 3, label.y = 30)
    
    #
    ####高通量数据处理####
    library(tidyverse)
    library(GEOquery)
    gset = getGEO('GSE167523', destdir=".", AnnotGPL = F, getGPL = F)
    gset[[1]]
    #通过pData函数获取分组信息
    pdata <- pData(gset[[1]])
    write.table(pdata, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
   exp0 <- read.table("exp_FPKM.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
   ids <- read.table("ids.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  
    exp <- as.data.frame(exp0)
   exp <- exp %>% mutate(ID=rownames(exp))
   ids <- ids %>% mutate(ID=rownames(ids))
   exp <- exp %>% inner_join(ids,by="ID") 
   exp <- exp[!duplicated(exp$Symbol),]#随机去重
   rownames(exp) <- exp$Symbol
   exp <- exp[,-(79:81)]#注意更换这个，删除最后两列
   write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
   
   boxplot(exp,outline=FALSE, notch=T)
   dev.off()    
   
   
   library(limma) 
   exp=normalizeBetweenArrays(exp)
   boxplot(exp,outline=FALSE, notch=T, las=2)
   range(exp)
   
   write.table(exp, file = "exp.normalize.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
   
   
   
   ####气泡矩阵图的绘制####
   pval2 <- read.csv("pval2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
   #注释：header=T表示数据中的第一行是列名，如果没有列名就用header=F
   data <- pval2
   data$group <- rownames(data)
   data$group = factor(data$group,
                           levels = c("HO",'SS','NASH','HCC'))
   library(reshape2)
   library(ggplot2)                          
   #注释：package使用之前需要调用
   data_melt<-melt (data,id.vars = "group")
   names(data_melt) = c('group', 'module', 'Value')
   ####注释：melt()函数把表格中的宽数据变成长数据####
   p<-ggplot(data_melt, aes(x = group, y = module, size = Value, color=group)) + geom_point()
   
   p
   #图片美化
   p<-ggplot(data_melt, aes(x = group, y = module, size = Value, color=group)) + geom_point()+
     theme(panel.background = element_blank(),
           panel.grid.major = element_line(colour = "gray"),
           panel.border = element_rect(colour="black",fill=NA))
p   
   
   
####多组GSEA分析####

####最新####
#limma差异分析
S1 <- read.table("S1group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
S2 <- read.table("S2group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exp <- read.table("normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exp <- exp[,rownames(group)]

library(limma) 



data<-exp
group_list<-read.table(file="group.txt",header=T,sep="\t",row.names=1)
#new_group<- group_list[order(group_list[,1]),]  #对分组数据进行排序，按照数据框的第一列升序排序
group<-new_group#取第二列，这样处理得到的就是每一个位置所对应的样本属性信息

group <- group_list
suppressMessages(library(limma))


design <- model.matrix(~0+factor(group$group))
design <- model.matrix(~0+factor(group))

rownames(design) <- rownames(group)

colnames(design)=c("group1",'group2')


#design数据表要符合都是数字形式0、1，data看情况要转置

contrast.matrix<-makeContrasts("group1-group2",levels=design)

##step1
fit <- lmFit(exp,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 

write.table(nrDEG, file = "group1-group2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(design, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####其他文献模型####

group_list<-read.table(file="HC-NAFL.txt",header=T,sep="\t",row.names=1)
group_list<-read.table(file="NAFL-NASH.txt",header=T,sep="\t",row.names=1)





exp<-read.table(file="normalize.txt",header=T,sep="\t",row.names=1)
gene <- group_list$gene
exp2 <- exp[gene,]
exp2 <- as.data.frame(t(exp2))

exp2$riskscore <- -0.190052713*exp2$CCL2+0.675427091*exp2$FABP4+-0.91795346*exp2$MMP7+-0.181266023*exp2$MT1X
+-0.777987595*exp2$MT1F+0.135941039*exp2$MT1G+1.336138959*exp2$MT1H+-0.42412614*exp2$MT1M+0.605956212*exp2$SPP1+-1.012626846

exp2$riskscore <- -0.166311656*exp2$CCL2+0.494621105*exp2$FABP4+1.434246184*exp2$MMP7
+2.35905055*exp2$MT1F+1.159066248*exp2$MT1G+-2.090099881*exp2$MT1H+-0.843261703*exp2$MT1M+0.442308223*exp2$SPP1+-29.40911129


group<-read.table(file="group.txt",header=T,sep="\t",row.names=1)
exp2 <- as.data.frame(t(exp2))
exp2 <- exp2[,rownames(group)]
exp2 <- exp2[,rownames(re1)]

exp2$group <- group$group1
exp2$group <- re1$group1

library(tidyverse)
library(rstatix)
library(ggpubr)
dat <- exp2
p1 <- ggboxplot(dat,  x = "group", y = "riskscore",
                color = "group", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme(axis.title.x=element_blank())+ #删除x轴坐标名称
  labs(title = "GSE167523")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1


write.table(exp2, file = "模型SS-NASH.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


library(ROCR)


re2 <- read.table("re2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re2 <- read.table("HC-SS.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re3 <- read.table("HO-SS.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

re1 <- read.table("SS-HASH.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re2 <- read.table("SS-HCC.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re3 <- read.table("NASH-HCC.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

re1 <- as.data.frame(re1)
re2 <- ifelse(re$group1=='2',1,0)

xy <- x
xy <- xy[rownames(re1),]
xy <- as.data.frame(xy)

re=cbind(re1 ,xy)
re$NASH <- re2
re <- re[re$diagnosis!='1',]
xy <- xy[xy$diagnosis!='1',]
re$PDZK1 <- xy$PDZK1
re$FADS2 <- xy$FADS2



pred_min1 <- prediction(re1[,2], re1[,1])
auc_min1 = performance(pred_min1,"auc")@y.values[[1]]

pred_min2 <- prediction(re2[,1], re2[,2])
auc_min2 = performance(pred_min2,"auc")@y.values[[1]]
pred_min3 <- prediction(re3[,2], re3[,1])
auc_min3 = performance(pred_min3,"auc")@y.values[[1]]

#求得AUC值
perf_min1 <- performance(pred_min1,"tpr","fpr")
perf_min2 <- performance(pred_min2,"tpr","fpr")
perf_min3 <- performance(pred_min3,"tpr","fpr")

plot(perf_min2,colorize=FALSE, col="red")
plot(perf_min1,colorize=FALSE, col="blue",add=T)
plot(perf_min3,colorize=FALSE, col="orange",add=T)

#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min1,3)))
text(0.8,0.3, labels = paste0("AUC = ",round(auc_min2,3)))
text(0.8,0.4, labels = paste0("AUC = ",round(auc_min3,3)))



legend("bottomright",
       c(paste0("AUC"),
         paste0("*SS vs NASH : ",round(auc_min2,3)), 
         paste0(" SS vs NASH  : ",round(auc_min1,3))),
       col=c('white',"red", "blue"),
       lty=1, lwd=2,bty = "n")  

# 加AUC值

pred_min1 <- prediction(re[,9], re[,3])
auc_min1 = performance(pred_min1,"auc")@y.values[[1]]

perf_min1 <- performance(pred_min1,"tpr","fpr")

plot(perf_min1,colorize=FALSE, col="blue")

lines(c(0,1),c(0,1),col = "gray", lty = 4 )

text(0.8,0.2, labels = paste0("AUC = ",round(auc_min1,3)))



####小鼠样本处理####
library(tidyverse)
library(GEOquery)
gset = getGEO('GSE40481', destdir=".", AnnotGPL = F, getGPL = F)
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
write.table(pdata, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

GPL=getGEO(filename = 'GSE40481_family.soft.gz')
# 提取信息（可以通用）
gpl=GPL@gpls[[1]]@dataTable@table
# 我只要ID和symbol
ids=gpl[,c(1,6)]
# 写出文件
write.table(ids,file = "ids.txt",sep = "\t",row.names=F,col.names = T)

head(ids)

#id转换
library(tidyverse)

exp <- exprs(gset[[1]])
ids$ID <- as.character(ids$ID)

exp <- as.data.frame(exp)
exp <- exp %>% mutate(ID=rownames(exp))
exp <- exp %>% inner_join(ids,by="ID") 
exp <- exp[!duplicated(exp$`Gene Symbol`),]#随机去重
rownames(exp) <- exp$`Gene Symbol`
exp <- exp[,-(171:172)]#注意更换这个，删除最后两列


#一个探针对应于多个基因的情况
GSE164760$symbol <- data.frame(sapply(GSE164760$symbol,
                                      function(x)unlist(strsplit(x,'///'))[1]),
                               stringsAsFactors = F)[,1]

exp <- aggregate(exp[,1:51],by=list(exp$ILMN_Gene),mean,na.rm= TRUE)#排除最后的一列group
rownames(exp) <- exp$Group.1
exp <- exp[,-1]



exp$gene <- rownames(exp)
exp <- exp[,-33]

exp$gene <- toupper(exp$gene)
rownames(exp) <- exp$gene

write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp, file = "exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



boxplot(exp,outline=FALSE, notch=T, las=2)
dev.off()
##数据矫正
library(limma) 
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
#判断是否log2取值

qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { exp[which(exp <= 0)] <- NaN
exp <- log2(exp)
print("log2 transform finished")}else{print("log2 transform not needed")}


range(exp)
dev.off()

write.table(exp, file = "exp.normalize_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


gene <- read.table("gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- gene$gene
data <- exp[gene,]
rownames(data) <- gene

write.table(data, file = "new.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


new <- read.table("new.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

new <- t(new)
new <- as.data.frame(new)
lasso.prob <- predict(cv_fit, newx=new , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
lasso.prob <- lasso.prob[rownames(group),]
re <- cbind(group,lasso.prob)

write.table(re, file = "re.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

group1 <- read.table("group1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group$source_name_ch1 <- group$characteristics_ch1.1
re1 <- re[rownames(group1),]
re1$source_name_ch1 <- group1$source_name_ch1

re2 <- re[rownames(group2),]
re2$source_name_ch1 <- group2$source_name_ch1



library(ggpubr) 
p1 <- ggboxplot(re2,  x = "characteristics_ch1.1", y = "s1",
                color = "characteristics_ch1.1", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  #删除x轴坐标名称
  labs(title = "GSE40481")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1



bxp <- ggboxplot(
  anxiety, x = "group", y = "score",
  color = "time", palette = "jco"
)
# 添加显著性标记
stat.test <- stat.test %>% add_xy_position(x = "group")
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)


p1 + stat_pvalue_manual(
  dat_p, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE, tip.length = 0
)


CYP7A1 <- exp["CYP7A1",,drop=F]
CYP7A1 <- t(CYP7A1)
CYP7A1 <- cbind(CYP7A1,group)
p1 <- ggboxplot(CYP7A1,  x = "source_name_ch1", y = "CYP7A1",
                color = "source_name_ch1", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  #删除x轴坐标名称
  labs(title = "GSE83596")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1

library(tidyverse)
library(rstatix)
library(ggpubr)


data <- read.table("re.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group1 <- read.table("group3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

data <- data[rownames(group1),]
data$type <- group1$type
colnames(data) <- c("1","2","risk score","s2","type")

p1 <- ggboxplot(data,  x = "type", y = "risk score",
                color = "type", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+
  theme(axis.title.x=element_blank())+#删除x轴坐标名称
  labs(title = "GSE83596")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))+
  theme(text = element_text(size = 20))
p1


stat.test <- data %>%
  group_by(type) %>%
  pairwise_t_test(
    s1 ~ type, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test

bxp <- ggboxplot(
  data, x = "group", y = "score",
  color = "time", palette = "jco"
)
# 添加显著性标记
stat.test <- stat.test %>% add_xy_position(x = "group")
p1 + stat_pvalue_manual(
  dat_p, label = "p.adj.signif", 
  step.increase = 0.08
)+stat_compare_means()


library(tidyverse)
library(rstatix)
library(ggpubr)
data("iris")
#按照Species分组，每组抽取一个样本
iris%>%
  sample_n_by(Species,size = 1)
#按照group分组，统计变量weight的标准差和方差
iris%>%
  group_by(Species)%>%
  get_summary_stats(Sepal.Length,
                    type = "mean_sd")
#使用ANOVA test比较多组的平均值
iris.res.aov<-iris%>%
  anova_test(Sepal.Length~Species)
iris.res.aov
#组件t检验
iris_test<-iris%>%
  pairwise_t_test(Sepal.Length~Species,
                  p.adjust.method = "bonferroni")
iris_test
#p值得坐标轴参数
iris_p<-iris_test%>%
  add_xy_position(x = "Species")
iris_p
#可视化
ggboxplot(iris,x = "Species",y = "Sepal.Length")+
  stat_pvalue_manual(iris_p,
                     label = "p.adj",
                     tip.length = 0,
                     step.increase = 0.1)+
  labs(
    subtitle = get_test_label(iris.res.aov,detailed = TRUE),
    caption = get_pwc_label(iris_p))



dat%>%
  sample_n_by(group,size = 1)
#按照group分组，统计变量weight的标准差和方差
dat%>%
  group_by(group)%>%
  get_summary_stats(colnames(dat[2]),
                    type = "mean_sd")
#使用ANOVA test比较多组的平均值
dat.res.aov<-data%>%
  anova_test(s1~type)
dat.res.aov
#组件t检验
dat_test<-data%>%
  pairwise_t_test(s1~type,
                  p.adjust.method = "bonferroni")
dat_test
#p值得坐标轴参数
dat_p<-dat_test%>%
  add_xy_position(x = "type")
dat_p
#可视化
ggboxplot(data,x = "type",y = colnames(data[3]),
          color = "type", palette = "Dark2",
          add = "jitter")+
  stat_pvalue_manual(dat_p,
                     label = "p.adj",
                     tip.length = 0,
                     step.increase = 0.1)+            #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  labs(
    subtitle = get_test_label(dat.res.aov,detailed = TRUE),
    caption = get_pwc_label(dat_p))


p1 <- ggboxplot(dat,  x = "group", y = colnames(dat[1]),
                color = "group", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme(axis.title.x=element_blank())+ #删除x轴坐标名称
  labs(title = "GSE167523")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1

####PDZK1在NASH中表达####
data0 <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

PDZK1 <- data0["PDZK1",]

PDZK1 <- as.data.frame(t(PDZK1))

PDZK1 <- PDZK1[rownames(group),,drop=F]

PDZK1 <- cbind(PDZK1,group)

data <- subset(PDZK1,group1=="HC"|group1=="NASH"|group1=="SS")

library(tidyverse)
library(rstatix)
library(ggpubr)
data("iris")
#按照Species分组，每组抽取一个样本
iris%>%
  sample_n_by(Species,size = 1)
#按照group分组，统计变量weight的标准差和方差
iris%>%
  group_by(Species)%>%
  get_summary_stats(Sepal.Length,
                    type = "mean_sd")
#使用ANOVA test比较多组的平均值
iris.res.aov<-iris%>%
  anova_test(Sepal.Length~Species)
iris.res.aov
#组件t检验
iris_test<-iris%>%
  pairwise_t_test(Sepal.Length~Species,
                  p.adjust.method = "bonferroni")
iris_test
#p值得坐标轴参数
iris_p<-iris_test%>%
  add_xy_position(x = "Species")
iris_p
#可视化
ggboxplot(iris,x = "Species",y = "Sepal.Length")+
  stat_pvalue_manual(iris_p,
                     label = "p.adj",
                     tip.length = 0,
                     step.increase = 0.1)+
  labs(
    subtitle = get_test_label(iris.res.aov,detailed = TRUE),
    caption = get_pwc_label(iris_p))

data <- GSVA
data <- as.data.frame(t(data))
group <- group0[rownames(data),]
data <- cbind(data,group)
col <- c("HALLMARK_MYC_TARGETS_V1","group1")
dat <- data[,col]
dat$group1 <- factor(dat$group1,levels = c("HC","HO","SS","NASH","HCC"))
colnames(dat) <- c("HALLMARK_MYC_TARGETS_V1","group")
dat%>%
  sample_n_by(group,size = 1)
#按照group分组，统计变量weight的标准差和方差
dat%>%
  group_by(group)%>%
  get_summary_stats(colnames(dat[1]),
                    type = "mean_sd")
#使用ANOVA test比较多组的平均值
dat.res.aov<-dat%>%
  anova_test(HALLMARK_MYC_TARGETS_V1~group)
dat.res.aov
#组件t检验
dat_test<-dat%>%
  pairwise_t_test(HALLMARK_MYC_TARGETS_V1~group,
                  p.adjust.method = "bonferroni")
dat_test
#p值得坐标轴参数
dat_p<-dat_test%>%
  add_xy_position(x = "group")
dat_p
#可视化
ggboxplot(dat,x = "group",y = "HALLMARK_MYC_TARGETS_V1",
          color = "group", palette = "Dark2",
          add = "jitter")+
  stat_pvalue_manual(dat_p,
                     label = "p.adj.signif",
                     tip.length = 0,
                     step.increase = 0.1)+            #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  labs(
    subtitle = get_test_label(dat.res.aov,detailed = T),
    caption = get_pwc_label(dat_p))


p1 <- ggboxplot(dat,  x = "group", y = "KEGG_FATTY_ACID_METABOLISM",
                color = "group", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme(axis.title.x=element_blank())+ #删除x轴坐标名称
  labs(title = "GSE167523")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1




GSE130970 <- read.table("GSE130970.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
GSE130970group <- read.table("GSE130790group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


PDZK1_2 <- GSE130970["PDZK1",]
PDZK1_2 <- as.data.frame(t(PDZK1_2))
PDZK1_2 <- PDZK1_2[rownames(GSE130970group),,drop=F]
data <- cbind(PDZK1_2,GSE130970group)
data$`fibrosis stage` <- factor(data$`fibrosis stage`,levels = c("0","1","2","3","4"))
data$`cytological ballooning grade` <- factor(data$`cytological ballooning grade`,levels = c("0","1","2"))
data$`nafld activity score` <- factor(data$`nafld activity score`,levels = c("0","1","2","3","4","5","6"))
data$`steatosis grade` <- factor(data$`steatosis grade`,levels = c("0","1","2","3"))

p1 <- ggboxplot(data,  x = "steatosis grade", y = "PDZK1",
                color = "steatosis grade", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme()+ #删除x轴坐标名称
  labs(title = "GSE130970")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1
p1 <- ggboxplot(data,  x = "cytological ballooning grade", y = "PDZK1",
                color = "cytological ballooning grade", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme()+ #删除x轴坐标名称
  labs(title = "GSE130970")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1



GSE135251 <- read.table("GSE135251.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
GSE135251group <- read.table("GSE135251group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

PDZK1 <- GSE135251["PDZK1",]
PDZK1 <- as.data.frame(t(PDZK1))
PDZK1 <- PDZK1[rownames(GSE135251group),,drop=F]
data <- cbind(PDZK1,GSE135251group)
data$`fibrosis stage:ch1` <- factor(data$`fibrosis stage:ch1`,levels = c("0","1","2","3","4"))
data$`group in paper:ch1` <- factor(data$`group in paper:ch1`,levels = c("control","NAFL","NASH_F0-F1","NASH_F2","NASH_F3","NASH_F4"))
data$`nas score:ch1` <- factor(data$`nas score:ch1`,levels = c("0","1","2","3","4","5","6","7","8"))

p1 <- ggboxplot(data,  x = "nas score:ch1", y = "PDZK1",
                color = "nas score:ch1", palette = "Set3",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme()+ #删除x轴坐标名称
  labs(title = "GSE135251")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1

GSE167523 <- read.table("GSE167523.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
GSE167523group <- read.table("GSE167523group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

PDZK1 <- GSE167523["PDZK1",]
PDZK1 <- as.data.frame(t(PDZK1))
PDZK1 <- PDZK1[rownames(GSE167523group),,drop=F]
data <- cbind(PDZK1,GSE167523group)
data$`disease subtype`
p1 <- ggboxplot(data,  x = "group1", y = "PDZK1",
                color = "group1", palette = "Set2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme()+ #删除x轴坐标名称
  labs(title = "")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1


#2024.7.3
#PDZK1和评分之间的相关性

data0 <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
GSVA <- read.table("gs.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("phe2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- na.omit(group2)

PDZK1 <- PDZK1[rownames(group2),]
data <- cbind(PDZK1,group2)

data <- subset(data,group1!="HC")

data <- subset(data,group1!="NASH")

PDZK1 <- data0["PDZK1",]

PDZK1 <- as.data.frame(t(PDZK1))

PDZK1 <- PDZK1[rownames(group),,drop=F]

PDZK1 <- cbind(PDZK1,group)

data <- subset(PDZK1,group1=="HC"|group1=="NASH"|group1=="SS")

sample <- rownames(data)
GSVA <- as.data.frame(t(GSVA))
gs.exp <- GSVA[sample,]
gs.exp <- na.omit(gs.exp)
sample <- rownames(gs.exp)
PDZK1 <- PDZK1[sample,]
data <- cbind(PDZK1,gs.exp)
data$KEGG_FATTY_ACID_METABOLISM
data$HALLMARK_MYC_TARGETS_V1

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))


b <- ggplot(data, aes(x =PDZK1 , y = KEGG_FATTY_ACID_METABOLISM))
b + geom_point()+
  geom_smooth( color="#2b83ba",method = "lm") +
  scale_color_manual(values = c("#2b83ba", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#E7B800", "#FC4E07"))+ 
  
  ggpubr::stat_cor()+
  labs(x="Relative PDZK1 Expression")

p1 <- ggboxplot(data,  x = "group1", y = "Age",
                color = "group1", palette = "Set2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme()+ #删除x轴坐标名称
  labs(title = "")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1

#2024.7.4
data <- subset(PDZK1,group1=="HC"|group1=="NASH")

data$type <- ifelse(data$group1=="HC","0","1")
library(ROCR)
pred_min <- prediction(data[,1], data[,4])
auc_min = performance(pred_min,"auc")@y.values[[1]]
#求得AUC值
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
# 加AUC值
library(pROC)
dfroc1<- roc(data[,4], data[,1])

plot(dfroc1,
     #col="black" ##曲线颜色
     #identity.col="blue"##对角线颜色
     print.auc=TRUE,print.auc.x=0.5,print.auc.y=0.5,  # 图像上输出AUC值,坐标为(x，y)
     auc.polygon=TRUE, auc.polygon.col="skyblue",  # 设置ROC曲线下填充及颜色
     max.auc.polygon=TRUE,   # 填充整个图像,ROC去曲线上填充
     ##partial.auc=c(1, 0.8), partial.auc.focus="sp",  # 范围内部分AUC；sp:横坐标；se：纵坐标
     grid=c(0.1, 0.2),grid.col=c("green", "red"), # 设置间距为0.1，0.2，线条颜色
     print.thres=F, print.thres.col="black",thresholds="best",  # 图像上输出最佳截断值,thresholds赋值给best
     smooth=TRUE,  #使ROC曲线变得光滑
     reuse.auc=F)# F不变，T是使横轴从0到1，表示为1-特异度

#2024.7.8
data0 <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gs.exp <- read.table("gs.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("phe2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group0 <- group0[colnames(gs.exp),]


data3 <- as.data.frame(t(data3))

data <- data[rownames(group),]
data <- cbind(data3,group3)

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)

p1 <- ggboxplot(data,  x = "group1", y = "KEGG_FATTY_ACID_METABOLISM",
                color = "group1", palette = "Set2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme()+ #删除x轴坐标名称
  labs(title = "")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p1

data <- subset(data,group1=="HC"|group1=="NASH"|group1=="SS")



#2024.7.9
data0 <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gs.exp <- read.table("gs.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gs.exp2 <- read.table("gs.exp2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
GSVA <- rbind(gs.exp,gs.exp2)
group2 <- read.table("phe2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
PDZK1 <- data0["PDZK1",]
PDZK1 <- as.data.frame(t(PDZK1))
group$PDZK1 <- PDZK1$PDZK1
group0 <- group
group1 <- subset(group0,group1=="SS"|group1=="NASH")
group1 <- subset(group0,group1=="NASH")
group3 <- subset(group0,group1=="NASH"|group1=="HCC")

group2 <- subset(group0,group1=="HCC")


group1$group <- ifelse(group1$PDZK1>=median(group1$PDZK1),"high","low")
group2$group <- ifelse(group2$PDZK1>=median(group2$PDZK1),"high","low")

data <- GSVA
data1 <- data0[,rownames(group1)]
data2 <- data[,rownames(group2)]
data3 <- data[,rownames(group3)]

library(limma)
data<-data2


#new_group<- group_list[order(group_list[,1]),]  #对分组数据进行排序，按照数据框的第一列升序排序
group <- group2#取第二列，这样处理得到的就是每一个位置所对应的样本属性信息



suppressMessages(library(limma))


design <- model.matrix(~0+factor(group$group))


rownames(design) <- rownames(group)

colnames(design)=c("high","low")


#design数据表要符合都是数字形式0、1，data看情况要转置

contrast.matrix<-makeContrasts("high-low",levels=design)

##step1
fit <- lmFit(data2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 

write.table(nrDEG, file = "HCC.GSVA_high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(group1, file = "NAFLDgroup.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(group2, file = "HCCgroup.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#火山图绘制

#火山图可视化

library(ggplot2)
library(ggrepel)
# 设置工作目录


# 整理数据集
# 参数'./resource/dataset.txt'，表示载入E:/R/WorkSpace/baimoc/visualization/resource/dataset_heatmap.txt
dataset <- read.table("NAFLD.GSVA_high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 设置pvalue和logFC的阈值
cut_off_adj.P.Val = 0.05
cut_off_logFC = 0.2
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
# 这里的change列用来设置火山图点的颜色
dataset$change = ifelse(dataset$adj.P.Val < cut_off_adj.P.Val & abs(dataset$logFC) >= cut_off_logFC, 
                        ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')

table(dataset$change)


# 绘制火山图
ggplot(
  #设置数据
  dataset, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#66C3A5", "#d2dae2","#FD8D62"))+
  
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_adj.P.Val),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (Adjust p-value)")+
  theme_bw()+
  
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

# 将需要标记的基因放置在label列
#添加基因标签信息
dataset$Label = ""   #新加一列label
dataset <- dataset[order(dataset$adj.P.Val), ]   #对差异基因的p值进行从小到大的排序
dataset$Gene <- rownames(dataset)
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(dataset$Gene[which(dataset$change == "Up")], 5)
#低表达的基因中，选择fdr值最小的5个


down.genes <- head(dataset$Gene[which(dataset$change == "Down")], 5)
#将up.genes和down.genes合并，并加入到Label中
dataset.top5.genes <- c(as.character(up.genes), as.character(down.genes))
dataset.top5.genes <- c("HALLMARK_MYC_TARGETS_V1","KEGG_FATTY_ACID_METABOLISM")
dataset$Label[match(dataset.top5.genes, dataset$Gene)] <- dataset.top5.genes

p4<-ggplot(
  #设置数据
  dataset, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c( "#d2dae2","#FD8D62"))+
  
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_adj.P.Val),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (Adjust p-value)")+
  theme_bw()+
  
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )+labs(title = "HCC vs HC")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))

p4

p1 <- p4+geom_text_repel(data = dataset, aes(x = dataset$logFC, 
                                            y = -log10(dataset$adj.P.Val), 
                                            label = Label),
                        size = 3,box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE)+labs(title = "High vs Low")+
  theme(plot.title = element_text(family = "serif"))

p1

#2024.7.10
#随机森林模拟
data <- as.data.frame(t(data2))
data$group <- group2$group
data$group <- factor(data$group)

library(randomForest)
library(ggplot2)
library(gridExtra)
# 构建随机森林模型
rf_model <- randomForest(group ~ ., data = data, ntree = 500, mtry = 2, importance = TRUE)
varImpPlot(rf_model, main = "Variable Importance")


# 提取变量重要性
importance_values <- importance(rf_model)

# 查看变量重要性
print(importance_values)

# 将变量重要性转换为数据框
importance_df <- as.data.frame(importance_values)

# 添加变量名列
importance_df$Variable <- rownames(importance_df)

# 查看转换后的数据框
print(importance_df)

# 确认变量重要性列名
colnames(importance_df)

if ("MeanDecreaseGini" %in% colnames(importance_df)) {
  importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]
  selected_metric <- "MeanDecreaseGini"
} else if ("MeanDecreaseAccuracy" %in% colnames(importance_df)) {
  importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  selected_metric <- "MeanDecreaseAccuracy"
} else if ("IncNodePurity" %in% colnames(importance_df)) {
  importance_df <- importance_df[order(importance_df$IncNodePurity, decreasing = TRUE), ]
  selected_metric <- "IncNodePurity"
} else {
  stop("No suitable metric found in importance output.")
}

# 查看排序后的变量重要性
print(importance_df)

# 使用 ggplot2 绘制变量重要性图
importance_plot <- ggplot(importance_df, aes(x = reorder(Variable, .data[[selected_metric]]), y = .data[[selected_metric]])) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Variable Importance", x = "Variables", y = selected_metric) +
  theme_minimal()

print(importance_plot)


# 筛选前n个最重要的指标
top_n <- 50
important_vars <- importance_df$Variable[1:top_n]

importance_df2 <- importance_df[1:top_n,]


importance_plot <- ggplot(importance_df2, aes(x = reorder(Variable, .data[[selected_metric]]), y = .data[[selected_metric]])) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Variable Importance", x = "Variables", y = selected_metric) +
  theme_minimal()

print(importance_plot)







write.table(importance_df, file = paste0("随机森林/","LIHC.fig/","ssGSEA_importance",".txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(importance_df, file = paste0("随机森林/","LIHC.fig/","GSVA_importance",".txt"),sep = "\t",row.names = T,col.names = NA,quote = F)


#KEGG和GOBP功能富集
library(tidyverse)
library(AnnotationDbi)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)



hh <- read.table("KEGG_NAFLD.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
hh <- hh[1:30,]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
ggplot(hh,aes(y=order,x=Count,fill=PValue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
  labs(
    x = "Gene numbers", 
    y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*PValue))+# 修改点的大小
  scale_color_gradient(low="#66C3A5",high = "#FD8D62")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Number",y="Pathways")+
  theme_bw()

hh <- read.table("GOBP.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
ggplot(hh,aes(y=order,x=Count,fill=PValue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
  labs(
    x = "Gene numbers", 
    y = "GO term")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

