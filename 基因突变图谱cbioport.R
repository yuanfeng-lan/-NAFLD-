####安装R包####
site = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = site)

current_packages = rownames(installed.packages())

needed_packages = c("maftools", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg19")

for (i in needed_packages) {
  if (!i %in% current_packages)
    BiocManager::install(i, update = F, site_repository=site)
  current_packages = rownames(installed.packages())
}

needed_packages2 <- c("mclust", "NMF", "barplot3d", "knitr")

for (j in needed_packages2) {
  if (!j %in% current_packages)
    install.packages(j)
  current_packages = rownames(installed.packages())
}

library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
library(barplot3d)
library(openxlsx)



#读入MAF文件
luad <- read.maf(maf="data_mutations.maf", clinicalData="group.txt")

luad <- read.maf(maf="data_mutations.maf")

luad # 直接输入MAF对象可以查看MAF文件的基本信息
getSampleSummary(luad) # 显示样品的统计
getGeneSummary(luad) # 显示基因的统计
getClinicalData(luad) # 显示样品关联的临床数据
getFields(luad) # 显示MAF文件中的所有字段
write.mafSummary(maf=luad, basename="luad") # 将详细统计结果输出到文件

#汇总统计表
plotmafSummary(maf=luad, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE,)

sample <- getSampleSummary(luad)

gene <- getGeneSummary(luad)
gene <- gene[1:20,]
gene <- gene$Hugo_Symbol

write.table(sample, file = "sample.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#瀑布图，top10表示展示前10个显著突变的基因
oncoplot(maf=luad, top=20,draw_titv = TRUE)
#添加标签以及标签颜色
#Color coding for FAB classification
fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)


#Color coding for FAB classification
fabcolors = RColorBrewer::brewer.pal(n = 3,name = 'Spectral')
names(fabcolors) = c("high", "low",'NA')
fabcolors = list(risk = fabcolors)

print(fabcolors)
#> $FAB_classification
#>        M0        M1        M2        M3        M4        M5        M6        M7 
#> "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#E6F598" "#ABDDA4" "#66C2A5" "#3288BD"

MYC_gene <- read.table("HALLMARK_MYC_TARGETS_V1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

aml_genes = c("TP53", "WT1", "PHF6", "DNMT3A", "DNMT3B", "TET1", "TET2", "IDH1", "IDH2", "FLT3", "KIT", "KRAS", "NRAS", "RUNX1", "CEBPA", "ASXL1", "EZH2", "KDM6A")

aml_genes = c('TP53','PIK3CA','LRP1B','KRAS','APC','FAT4','KMT2D','KMT2C','BRAF')
aml_genes <- MYC_gene$gene
oncoplot(maf = luad, clinicalFeatures = 'risk',)


oncoplot(
  maf = luad, top=20,
  clinicalFeatures = 'risk',
  sortByAnnotation = TRUE,
  annotationColor = fabcolors
)

#分组瀑布图
# 从临床数据中提取性别对应的"Tumor_Sample_Barcode"
clin <- read.table("PDZK1.txt", header=T, sep="\t")

clin.low <- subset(clin, PDZK1=="low")$Tumor_Sample_Barcode
clin.high <- subset(clin, risk=="high")$Tumor_Sample_Barcode
# 使用subsetMaf构建男性和女性的MAF对象
luad.low <- subsetMaf(maf=luad, tsb=clin.low, isTCGA=F)

luad.high <- subsetMaf(maf=luad, tsb=clin.high, isTCGA=F)

# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=luad.low, m2=luad.high, m1Name="low", m2Name="high", minMut=5)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="low_vs_high.tsv", quote=FALSE, row.names=FALSE, sep="\t")
results <- fvsm$results
genes <- results[1:20,]

genes <- genes$Hugo_Symbol

genes <- c("PCDH11Y", "HTATSF1", "WWC3", "DMD", "PLXNB3", "AMER", "MCF2", "FLNA")
coOncoplot(m1=luad.low, m2=luad.high, m1Name=" risk low", m2Name="risk high",genes = genes)

gene_results <- as.data.frame(gene_results)

gene_results <- results
rownames(gene_results) <- gene_results$Hugo_Symbol
gene_results <- gene_results[gene,]
write.table(gene_results, file = "gene_results.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#计算tmb值
tmb_table_wt_log = tmb(maf = luad)
write.table(tmb_table_wt_log, file = "TMB.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#


####TMB相关性分析####
library(datasets)
library(ggplot2)
library(ggpubr)
library(ggpmisc)


TMB <- read.table("exp.normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- read.table("gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- gene$gene
TMB <- TMB[gene,]
TMB <- t(TMB)
TMB <- as.data.frame(TMB)

my_data <- TMB

cor.test(my_data$total_perMB_log,my_data$PDZK1,method="pearson",use="complete.obs")


ggscatter(my_data, x = "NOL3", y = "FGF14",
          add = "reg.line", conf.int = T,    
          add.params = list(fill = "lightgray")
          
)+
  stat_cor(method = "pearson", 
           label.x = 10, label.y = 12)

#x\y两个参数调整标签的位置

####分组比较t检验####


library(ggpubr)

mydata <- read.table("TMB2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

p <- ggboxplot(mydata, x="PDZK1", y = "TMB", color = "PDZK1", palette = "jco", add = "jitter",  short.panel.labs = FALSE)

# 添加p值
p + stat_compare_means(method = "t.test",label.y=20)
# 显示p值但不显示方法
p + stat_compare_means(aes(label = ..p.format..),method = "t.test",label.x = 1.5)

# 只显示显著性水平
p + stat_compare_means(aes(label = ..p.signif..),method = "t.test",label.x = 1.5)

####卡方检验检测百分比或非连续变量组之间的差异####

Usage
chisq.test(x, y = NULL, correct = TRUE,
           p = rep(1/length(x), length(x)), rescale.p = FALSE,
           simulate.p.value = FALSE, B = 2000)



####生存分析####
setwd("survival")
surv <- exp_sur
surv <- read.table("survival.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
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
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
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
