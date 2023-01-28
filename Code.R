######################################################################
setwd("C:\\Users\\fancy\\Desktop\\TCGA")

rm(list = ls())
gc()

clinical <- read.csv("clinical.csv",header = T,row.names = 1)
rownames(clinical) <- gsub("TCGA","T",rownames(clinical))
rownames(clinical) <- gsub("-",".",rownames(clinical))

clinical <- clinical[!is.na(clinical$A1_OS),]

######################################################################
###TCGA STAD FPKM
exp <- read.table("TCGA_STAD_FPKM.txt",sep = "\t",header = T,row.names = 1)

library(dplyr)
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
TPM <- apply(exp,2,FPKM2TPM)
TPM <- data.frame(TPM)
colSums(TPM) 

for (i in 1:length(TPM)) {
  m <- colnames(TPM)[i]
  num <- as.numeric(substr(m,14,15))
  if (num %in% 1:9) {colnames(TPM)[i] = gsub(substr(m,1,4),"T",m)}
  else {colnames(TPM)[i] = gsub(substr(m,1,4),"N",m)}
}

#a <- substr(colnames(TPM),1,9)
#table(a[duplicated(a)])
colnames(TPM) <- substr(colnames(TPM),1,9)
log2TPM <- log2(TPM+1) 
test <- log2TPM[1:10,1:10]

save(clinical,log2TPM,file = "STAD_TPM_Clinical.Rdata")

######################################################################
boxplot(log2TPM[,1:20])

#install.packages("factoextra")
library(factoextra)

data <- t(log2TPM)
res.pca <- prcomp(data)
#fviz_screeplot(res.pca, addlabels = TRUE) 
condition = factor(substr(rownames(data),1,1))
fviz_pca_ind(res.pca, label="none", habillage=condition,
             addEllipses=TRUE, ellipse.level=0.95,
             palette = c("#999999", "#E69F00", "#56B4E9"))

library(pheatmap)         

data <- data.frame(t(log2TPM %>% sample_n(100)))
data$cluster <- substr(rownames(data),1,1)
data <- data[order(data$cluster),]
Type <- data.frame(data$cluster)
rownames(Type) <- rownames(data)
data <- t(data[,1:(ncol(data)-1)])

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(Type$data.cluster)))]
names(CluCol)=levels(factor(Type$data.cluster))
ann_colors[["cluster"]]=CluCol

pheatmap <- pheatmap(data,
                     annotation=Type,
                     annotation_colors = ann_colors,
                     color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
                     cluster_cols =F,
                     cluster_rows =F,
                     scale="row",
                     show_colnames=F,
                     show_rownames=F,
                     fontsize=6,
                     fontsize_row=5,
                     fontsize_col=6)

####################################################################

rm(list = ls())
gc()

load("STAD_TPM_Clinical.Rdata")

library(dplyr)
library(survival)

FATTY.genes <- read.delim("FATTY.ACID.METABOLISM.genes.txt",header=F)
TPM_T <- log2TPM %>% select(starts_with("T"))

merge <- merge(t(TPM_T),clinical,by="row.names")
merge <- merge %>% filter(A1_OS > 30)
rownames(merge) <- merge$Row.names
merge <- merge[,-1]
test <- merge[1:10,1:10]

TPM_T <- data.frame(t(merge[,1:length(rownames(TPM_T))]))
#save(clinical,log2TPM,TPM_T,file = "STAD_TPM_Clinical.Rdata")

FATTY.exp <- subset(TPM_T, rownames(TPM_T) %in% FATTY.genes$V1)

write.table(rbind(id=colnames(FATTY.exp),FATTY.exp),
            file="fatty243.exp.txt", sep="\t",quote=F,col.names = F)

###################################################################

rm(list = ls())
gc()

library(dplyr)
library(survival)

FATTY.exp <- read.table("fatty243.exp.txt",header=T,sep="\t",check.names=F,row.names=1)

coxcli <- read.csv("clinical.csv",header = T,row.names = 1)
rownames(coxcli) <- gsub("TCGA","T",rownames(coxcli))
rownames(coxcli) <- gsub("-",".",rownames(coxcli))

merge <- merge(t(FATTY.exp),coxcli,by="row.names")
rownames(merge) <- merge$Row.names
merge <- merge[,-1]
test <- merge[1:10,1:10]

################################################################################
###uni
merge$A1_OS <- merge$A1_OS/365
merge$A1_OS <- as.numeric(merge$A1_OS)
merge$A2_Event <- ifelse(merge$A2_Event == "Alive",c(0), c(1))

recu = c("A1_OS","A2_Event")
outTab=data.frame()
for(i in rownames(FATTY.exp)){
  cox <- coxph(Surv(A1_OS,A2_Event) ~ merge[,i], data = merge)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  if(coxP<0.05){
    recu=c(recu,i)
  }
}

write.csv(outTab,file="uniCox.csv")

uniSigExp=merge[,recu]
write.table(cbind(id=row.names(uniSigExp),uniSigExp),
            file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

################################################################################
###mul
uniSigExp <- read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
rt <- na.omit(uniSigExp)

multiCox=coxph(Surv(A1_OS,A2_Event) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
#save(multiCox,file = "multiCox.Rdata")

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

write.csv(outTab,file="multiCox.csv")

riskScore=predict(multiCox,type="risk",newdata=uniSigExp)

tab <- cbind(uniSigExp[,c("A1_OS","A2_Event",rownames(outTab))],riskScore)

write.table(rbind(ID=colnames(tab), tab),
            file="fatty.score.txt",
            sep="\t", quote=F, col.names=F)

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)

cluster <- cluster[order(cluster$riskScore),]
median <- median(cluster$riskScore)
cluster$cluster[cluster$riskScore <= median] = "low"
cluster$cluster[cluster$riskScore > median] = "high"

clu=data.frame(cluster[,"cluster"])
rownames(clu) <- rownames(cluster)
colnames(clu) <- "cluster"

library(pheatmap)         

data=cluster[,c(3:8,10)]
Type <- data.frame(data$cluster)
rownames(Type) <- rownames(data)
data <- t(data[,1:6])

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(Type$data.cluster)))]
names(CluCol)=levels(factor(Type$data.cluster))
ann_colors[["cluster"]]=CluCol

pheatmap<-pheatmap(data,
                   annotation=Type,
                   annotation_colors = ann_colors,
                   color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
                   cluster_cols =F,
                   cluster_rows =F,
                   scale="row",
                   show_colnames=F,
                   show_rownames=T,
                   fontsize=9,
                   fontsize_row=9,
                   fontsize_col=9)

################################################################################
#install.packages("survivalROC")
library(survivalROC)
library(survminer)
library(survival)

rt=read.table("fatty.score.txt",header=T,sep="\t",check.names=F,row.names=1)

sur_outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  
  m <- median(rt[,i])
  rt$cluster[rt[,i] < m] = "low"
  rt$cluster[rt[,i] > m] = "high"
  
  length=length(levels(factor(rt$cluster)))
  diff=survdiff(Surv(A1_OS,A2_Event) ~ cluster, data = rt)
  pValue=1-pchisq(diff$chisq, df=length-1)
  
  sur_outTab=rbind(sur_outTab,
                   cbind(id=i,pvalue=pValue)
  )
}

write.csv(sur_outTab,file="sur_outTab.csv")

if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

###i in colnames(rt[,3:ncol(rt)])
m <- median(rt[,"riskScore"])
rt$cluster[rt[,"riskScore"] < m] = "low"
rt$cluster[rt[,"riskScore"] > m] = "high"

fit <- survfit(Surv(A1_OS,A2_Event) ~ cluster, data = rt)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
ggsurvplot(fit,
           data=rt,
           conf.int=F,
           pval=pValue,
           pval.size=6,
           legend.title="Cluster",
           legend.labs=levels(factor(rt[,"cluster"])),
           legend = c(0.8, 0.8),
           font.legend=10,
           xlab="Time(years)",
           break.time.by = 1,
           palette = bioCol,
           surv.median.line = "hv",
           risk.table=T,
           cumevents=F,
           risk.table.height=.25)
#########################################################################

rm(list = ls())
gc()

library(limma)
library(dplyr)

load("STAD_TPM_Clinical.Rdata")

TPM_T <- t(TPM_T)

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)
cluster <- cluster[,8:9]

merge <- merge(cluster,TPM_T,by="row.names")
rownames(merge) <- merge[,1]
merge <- merge[,-1]
merge <- merge[order(merge$riskScore),-1]
test <- merge[1:10,1:10]

merge$cluster[merge$riskScore < 0.777] = "low"
merge$cluster[merge$riskScore > 1.288] = "high"
merge <- merge %>% filter(cluster == "high" | cluster == "low")
group <-  dplyr::select(merge,-"cluster")

group <- data.frame(t(group[,-1]))
test <- group[1:10,1:10]

list <- c(rep("low", 111), rep("high",111)) %>% factor(., levels = c("high", "low"), ordered = F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("high", "low")
df.fit <- lmFit(group, list) 

df.matrix <- makeContrasts(high - low, levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "BH")

diffsig = na.omit(tempOutput) 

foldChange = 1
padj = 0.05

All_diffSig <- diffsig[(diffsig$adj.P.Val < padj & 
                          (diffsig$logFC>foldChange | 
                             diffsig$logFC < (-foldChange))),]

write.table(rbind(id=colnames(All_diffSig),All_diffSig), 
            "h_l_diffsig.txt",sep="\t", quote=F,col.names=F)  

#####################################################################
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

rt=read.table("h_l_diffsig.txt",header=T,sep="\t",check.names=F)

#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#定义显示Term数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#气泡图
bub=dotplot(kk, 
            showCategory=showNum, 
            orderBy="GeneRatio", 
            split="ONTOLOGY", 
            color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)

#####################################################################
#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=20
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#气泡图
dotplot(kk, 
        showCategory=showNum, 
        orderBy="GeneRatio", 
        color=colorSel)

#####################################################################
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

load("STAD_TPM_Clinical.Rdata")

TPM_T <- t(TPM_T)

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)
cluster <- cluster[,8:9]

merge <- merge(cluster,TPM_T,by="row.names")
rownames(merge) <- merge[,1]
merge <- merge[,-1]

merge <- merge[order(merge$riskScore),-1]

merge$cluster[merge$riskScore < 0.777] = "low"
merge$cluster[merge$riskScore > 1.288] = "high"
merge <- merge %>% filter(cluster == "high" | cluster == "low")
merge <-  dplyr::select(merge,-"riskScore")

group <-  dplyr::select(merge,-"cluster")
group <- data.frame(t(group))
test <- group[1:10,1:10]

group <- data.matrix(group)
c2gmt <- getGmt("c2.cp.v7.5.1.symbols.gmt")
gene.set <- c2gmt[grep("^REACTOME", names(c2gmt)),]

gsvaResult=gsva(group, 
                gene.set, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

cluster=data.frame(merge[,"cluster"])
rownames(cluster) <- rownames(merge)
colnames(cluster) <- "cluster"

gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
test <- gsvaResult[1:10,1:10]

adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$cluster)

treat=gsvaCluster[gsvaCluster$cluster=="high",]
con=gsvaCluster[gsvaCluster$cluster=="low",]
data=rbind(con, treat)

Type=as.vector(data$cluster)
data=t(data[,-ncol(data)])
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
fit=lmFit(data, design)
cont.matrix=makeContrasts(high - low, levels=design)
fit2=contrasts.fit(fit, cont.matrix)
fit2=eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, "allDiffOut.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, "diffSigOut.txt", sep="\t", quote=F, col.names=F)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(allType)))]
names(CluCol)=levels(factor(allType))
ann_colors[["cluster"]]=CluCol

termNum=20
diffTermName=as.vector(rownames(diffSig))
diffLength=length(diffTermName)
if(diffLength<termNum){termNum=diffLength}
hmGene=diffTermName[1:termNum]
hmExp=data[hmGene,]

ann=data.frame(Type)
rownames(ann) <- rownames(t(data))
colnames(ann) <- "cluster"

pheatmap(hmExp, 
         annotation=ann,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         cluster_cols =F,
         show_colnames = F,
         gaps_col=as.vector(cumsum(table(Type))),
         scale="row",
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=10)

####################################################################
###ssGSEA

load("TCGA_STAD_sig_go_kegg.Rdata")

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)

cluster <- cluster[order(cluster$riskScore),]
cluster$cluster[cluster$riskScore < 0.777] = "low"
cluster$cluster[cluster$riskScore > 1.288] = "high"
merge <- cluster %>% filter(cluster == "high" | cluster == "low")

cluster=data.frame(merge[,"cluster"])
rownames(cluster) <- rownames(merge)
colnames(cluster) <- "cluster"

gsvaResult <- sig_go_kegg[,1:52]
rownames(gsvaResult) <- gsvaResult$ID

sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
test <- gsvaResult[1:10,1:10]

adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$cluster)

treat=gsvaCluster[gsvaCluster$cluster=="high",]
con=gsvaCluster[gsvaCluster$cluster=="low",]
data=rbind(con, treat)

Type=as.vector(data$cluster)
data=t(data[,3:52])
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
fit=lmFit(data, design)
cont.matrix=makeContrasts(high - low, levels=design)
fit2=contrasts.fit(fit, cont.matrix)
fit2=eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, "hall.allDiffOut.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, "hall.diffSigOut.txt", sep="\t", quote=F, col.names=F)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(allType)))]
names(CluCol)=levels(factor(allType))
ann_colors[["cluster"]]=CluCol

termNum=20
diffTermName=as.vector(rownames(diffSig))
diffLength=length(diffTermName)
if(diffLength<termNum){termNum=diffLength}
hmGene=diffTermName[1:termNum]
hmExp=data[hmGene,]

ann=data.frame(Type)
rownames(ann) <- rownames(t(data))
colnames(ann) <- "cluster"

pheatmap(hmExp, 
         annotation=ann,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         cluster_cols =F,
         show_colnames = F,
         gaps_col=as.vector(cumsum(table(Type))),
         scale="row",
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=10)

######################################################################
###IOBR
library(IOBR) 
library(dplyr)

load("STAD_TPM_Clinical.Rdata")

eset_stad <- TPM_T

#TCGA, RNAseq, log2(TPM+1)
sig_tme<-calculate_sig_score(
  pdata=NULL,
  eset=eset_stad,
  signature=signature_collection,
  method="pca", 
  adjust_eset = T,
  mini_gene_count = 2)

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)
cluster <- cluster[,8:9]

cluster$cluster[cluster$riskScore < 0.777] = "low"
cluster$cluster[cluster$riskScore > 1.288] = "high"
cluster <- cluster %>% filter(cluster == "high" | cluster == "low")

data <- sig_tme
rownames(data) <- data$ID
data <- merge(cluster,data,by="row.names")
rownames(data) <- data$Row.names

Diff <- c()
Diff_all <-c()
for(sub in colnames(sig_tme[,3:ncol(sig_tme)])){ 
  Data <- data[c("cluster",sub)]
  normoxic <- filter(Data,cluster=="low") %>% na.omit()
  normoxic[,2] = as.numeric(normoxic[,2])
  hypoxic <- filter(Data,cluster=="high") %>% na.omit()
  hypoxic[,2] = as.numeric(hypoxic[,2])
  
  L_shapiro <- c()
  H_shapiro <- c()
  Bartlett  <- c()
  try(L_shapiro <- shapiro.test(normoxic)$p.value,silent = T) 
  try(H_shapiro <- shapiro.test(hypoxic)$p.value,silent = T)
  Bartlett=bartlett.test(Data[,2]~factor(Data$cluster), data=Data)$p.value 
  
  if(sum(L_shapiro,na.rm = T) > 0.05 & sum(H_shapiro,na.rm = T) > 0.05 & sum(Bartlett,na.rm = T) > 0.05){
    diffmean <- mean(hypoxic[,2],na.rm = T)-mean(normoxic[,2],na.rm = T)
    Logmean <- mean(hypoxic[,2],na.rm = T)/mean(normoxic[,2],na.rm = T)
    Pvalue <- unlist(t.test(hypoxic[,2],normoxic[,2])$p.value)
    Method <- "T-test"
    Diff <- data.frame(diffmean = diffmean,Log2FC=log2(Logmean),Pvalue=Pvalue,Method=Method,row.names=sub)
  }else{
    diffmean <- mean(hypoxic[,2],na.rm = T)-mean(normoxic[,2],na.rm = T)
    Logmean <- mean(hypoxic[,2],na.rm = T)/mean(normoxic[,2],na.rm = T)
    Pvalue <- wilcox.test(hypoxic[,2],normoxic[,2])$p.value
    Method <- "Wilcox"
    Diff <- data.frame(diffmean = diffmean,Log2FC=log2(Logmean),Pvalue=Pvalue,Method=Method,row.names=sub)
  }
  Diff_all <- rbind(Diff_all,Diff)
}

Diff_all$adj.P.Val <- p.adjust(Diff_all$Pvalue,method="fdr")

write.csv(Diff_all,"hl_STAD_sig_collection_result.csv")

######################################################################
library("IOBR")
library("tidyHeatmap")
library(dplyr)
#install.packages("ggcorrplot")
library(ggcorrplot)

group <- data[,c("ID","cluster")]

res<-iobr_cor_plot(pdata_group           = group, 
                   id1                   = "ID",
                   feature_data          = feature.data, 
                   id2                   = "ID",
                   target                = NULL, 
                   group                 = "cluster", 
                   is_target_continuous  = FALSE, 
                   padj_cutoff           = 1, 
                   index                 = 1, 
                   category              = "signature", 
                   signature_group       = sig_group,
                   palette_box           = "paired1", 
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2, 
                   feature_limit         = 26, 
                   character_limit       = 30, 
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE, 
                   show_plot             = F, 
                   path                  = "Relevant-signatures")  

write.csv(res,"HL_relevant_signatures.csv")
#####################################################################
library("IOBR")

multicox <- t(eset_stad)[,c("ADH4","AKR1B1","CYP4A11","NEU2","SMPD3","ST6GALNAC3")]
multicox <- data.frame(multicox)
multicox$ID <- row.names(multicox)

res<-iobr_cor_plot(pdata_group           = multicox, 
                   id1                   = "ID", 
                   feature_data          = sig_tme, 
                   id2                   = "ID", 
                   target                = "ST6GALNAC3", 
                   group                 = "group3", 
                   is_target_continuous  = TRUE, 
                   padj_cutoff           = 1, 
                   category              = "signature", 
                   signature_group       = sig_group[c(2,4,14,31)], 
                   ProjectID             = "TCGA", 
                   palette_box           = "set2", 
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 26,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = F,
                   path                  = "ST6GALNAC3-relevant-signatures")

##############################################################################
rm(list = ls())
gc()

library(dplyr)
library(survival)

FATTY.exp <- read.table("fatty243.exp.txt",header=T,sep="\t",check.names=F,row.names=1)

coxcli <- read.csv("clinical.csv",header = T,row.names = 1)
rownames(coxcli) <- gsub("TCGA","T",rownames(coxcli))
rownames(coxcli) <- gsub("-",".",rownames(coxcli))

merge <- merge(t(FATTY.exp),coxcli,by="row.names")
rownames(merge) <- merge$Row.names
merge <- merge[,-1]
test <- merge[1:10,1:10]

merge$A17_Age <- ifelse(merge$A17_Age > 65,c(">65"), c("<=65"))
merge$A18_Sex <- ifelse(merge$A18_Sex == "MALE",c(1), c(0))
merge$A1_OS <- merge$A1_OS/365
merge$A1_OS <- as.numeric(merge$A1_OS)
merge$A2_Event <- ifelse(merge$A2_Event == "Alive",c(0), c(1))

recu = c("A1_OS","A2_Event")
outTab=data.frame()
for(i in colnames(merge)){
  cox <- coxph(Surv(A1_OS,A2_Event) ~ merge[,i], data = merge)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  if(coxP<0.05){
    recu=c(recu,i)
  }
}

write.csv(outTab,"unicox.csv")

uniSigExp=merge[,recu]
write.table(cbind(id=row.names(uniSigExp),uniSigExp),
            file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

################################################################################
rt <- na.omit(uniSigExp)
rt <- rt[,1:21]

multiCox=coxph(Surv(A1_OS,A2_Event) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

write.csv(outTab,file="multiCox.csv")

riskScore=predict(multiCox,type="risk",newdata=uniSigExp)

tab <- cbind(uniSigExp[,c("A1_OS","A2_Event","A17_Age","A3_T","A4_N",
                          rownames(outTab))],riskScore)

################################################################################
#nomogram calibration ROC C-index DCA NRI IDI
###nomogram
library(regplot)
library(survival)
library(timeROC)

rt <- tab
multiCox=coxph(Surv(A1_OS,A2_Event) ~ riskScore+A4_N+A3_T+A17_Age, data = rt)
multiCoxSum=summary(multiCox)

medianTrainRisk = median(rt$riskScore)
rt$risk <- ifelse(rt$riskScore > medianTrainRisk,"high","low")

rt$nomo_riskScore=predict(multiCox,type="risk",newdata=rt)

nomo_medianTrainRisk = median(rt$nomo_riskScore)
rt$nomo_risk <- ifelse(rt$nomo_riskScore > nomo_medianTrainRisk,"high","low")

write.table(rbind(ID=colnames(rt), rt),
            file="nomo.score.txt",
            sep="\t", quote=F, col.names=F)

#####################################################################
rt <- read.table("nomo.score.txt",header=T,sep="\t",check.names=F,row.names=1)       

###time AUC
mayo1 <- rt
mayo1$A2_Event <- ifelse(mayo1$A2_Event == 1,c(2), c(1))
ROC <- timeROC(T=mayo1$A1_OS,
               delta=mayo1$A2_Event, 
               marker=mayo1$nomo_riskScore,
               cause=2,weighting="marginal",
               times=quantile(mayo1$A1_OS,probs=seq(0.2,1,0.2)),
               iid=TRUE)
plotAUCcurve(ROC,conf.int=TRUE,col="red")

################################################################################
###nomogram
pbccox <-  coxph(formula = Surv(A1_OS,A2_Event) ~ A4_N+A3_T+A17_Age+riskScore,data=rt)
a1<-regplot(pbccox,observation=rt[3,],points=T,failtime=c(5,3,1),
            droplines = T,prfail = F,center=F,plots=c("boxes","boxes"))
#"no plot", "boxes", "bars" or "spikes"
#"no plot", "density", "boxes", "spikes", "ecdf", "bars", "boxplot", "violin" or "bean".

###calibration
library(foreign)
library(survival)
library(rms)
library(tidyr)
library(dplyr)
library(plyr)
#####################DFS#######################
cox1 <- cph(Surv(A1_OS,A2_Event) ~ riskScore+A4_N+A3_T+A17_Age,surv=T,x=T, y=T,time.inc =1,data=rt)
cal1 <- calibrate(cox1, cmethod="KM", method="boot", u=1, m=100, B=1000)
plot(cal1)
cox2 <- cph(Surv(A1_OS,A2_Event) ~ riskScore+A4_N+A3_T+A17_Age,surv=T,x=T, y=T,time.inc =3,data=rt)
cal2 <- calibrate(cox2, cmethod="KM", method="boot", u=3, m=100, B=1000)
plot(cal2)
cox3 <- cph(Surv(A1_OS,A2_Event) ~ riskScore+A4_N+A3_T+A17_Age,surv=T,x=T, y=T,time.inc =5,data=rt)
cal3 <- calibrate(cox3, cmethod="KM", method="boot", u=5, m=100, B=1000)
plot(cal3)

cal <- c(cal1,cal2,cal3)
plot(cal1,xlim = c(0,1),ylim = c(0,1),lwd=4,lty=1,pch = 19,
     xlab ="Nomogram-Predicted OS(%)",ylab="Observed OS(%)",conf.int = F,
     col = "#FF0000",
     sub=F)
lines(cal1[,c('mean.predicted','KM')],type = 'b',lwd = 4,col = "#FF0000",pch = 16)
plot(cal2,xlim = c(0,1),ylim = c(0,1),lwd=4,lty=3,pch = 19,
     xlab ="Nomogram-Predicted OS(%)",ylab="Observed OS(%)",conf.int = F,
     col = "#0000EE",add = T,
     sub=F)
lines(cal2[,c('mean.predicted','KM')],type = 'b',lwd = 4,col = "#0000EE",pch = 16)
plot(cal3,xlim = c(0,1),ylim = c(0,1),lwd=4,lty=3,pch = 19,
     xlab ="Nomogram-Predicted OS(%)",ylab="Observed OS(%)",conf.int = F,
     col = "#00EE00",add = T,
     sub=F)
lines(cal3[,c('mean.predicted','KM')],type = 'b',lwd = 4,col = "#00EE00",pch = 16)

legend("bottomright",c("1 year Survival","3 year Survival","5 year Survival"),
       col=c("#FF0000","#0000EE","#00EE00"),lty=1,lwd=2)

################################################################################
###ROC1
#install.packages("survivalROC")
library(survivalROC)
library(survminer)
rt=read.table("nomo.score.txt",header=T,sep="\t",check.names=F,row.names=1)

diff=survdiff(Surv(A1_OS,A2_Event) ~ rt$nomo_risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

fit <- survfit(Surv(A1_OS,A2_Event) ~ nomo_risk, data = rt)

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=T,
                   pval=pValue,
                   pval.size=6,
                   legend.title="Risk",
                   legend.labs=c("0", "1"),
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette=c("red", "blue"),
                   risk.table=TRUE,
                   risk.table.title="",
                   risk.table.col = "strata",
                   risk.table.height=.25)

print(surPlot)
##################DFS#############
mayo1=read.table("nomo.score.txt",header=T,sep="\t",check.names=F,row.names=1)
cutoff <- 5
######model 1###########
Mayo4.1= survivalROC(Stime=mayo1$A1_OS,  
                     status=mayo1$A2_Event,      
                     marker = mayo1$A17_Age,     
                     predict.time =  cutoff, method="KM")
######model 2#############
Mayo4.2= survivalROC(Stime=mayo1$A1_OS,  
                     status=mayo1$A2_Event,      
                     marker = mayo1$A3_T,     
                     predict.time = cutoff, method="KM")
######model 3#############
Mayo4.3= survivalROC(Stime=mayo1$A1_OS,  
                     status=mayo1$A2_Event,      
                     marker = mayo1$A4_N,     
                     predict.time = cutoff, method="KM")
########nomogram
Mayo4.4= survivalROC(Stime=mayo1$A1_OS,  
                     status=mayo1$A2_Event,      
                     marker = mayo1$riskScore,     
                     predict.time = cutoff, method="KM")
########nomogram
Mayo4.5= survivalROC(Stime=mayo1$A1_OS,  
                     status=mayo1$A2_Event,      
                     marker = mayo1$nomo_riskScore,     
                     predict.time = cutoff, method="KM")

plot(Mayo4.5$FP, Mayo4.5$TP, ## x=FP,y=TP
     type="l",col="#0000EE",lwd=4, 
     xlim=c(0,1), ylim=c(0,1),   
     xlab=("FP"), 
     ylab="TP",
     main="Predicted survival")
abline(0,1,col="black",lty=2)
## nomogram
lines(Mayo4.4$FP, Mayo4.4$TP, type="l",col="#EE0000",lwd=4,xlim=c(0,1), ylim=c(0,1))
lines(Mayo4.3$FP, Mayo4.3$TP, type="l",col="#FF00FF",lwd=4,xlim=c(0,1), ylim=c(0,1))
lines(Mayo4.2$FP, Mayo4.2$TP, type="l",col="#FF0000",lwd=4,xlim=c(0,1), ylim=c(0,1))
lines(Mayo4.1$FP, Mayo4.1$TP, type="l",col="#00FF00",lwd=4,xlim=c(0,1), ylim=c(0,1))
legend(0.7,0.3,c(paste("Nomogram AUC=",round(Mayo4.5$AUC,3)),
                 paste("LMAGs AUC=",round(Mayo4.4$AUC,3)),
                 paste("N AUC=",round(Mayo4.3$AUC,3)),
                 paste("T AUC=",round(Mayo4.2$AUC,3)),
                 paste("age AUC=",round(Mayo4.1$AUC,3))
),
x.intersp=1, y.intersp=1,
lty= 1 ,lwd= 2,col=c("#0000EE","#EE0000","#FF00FF","#FF0000","#00FF00"),
bty = "n",
seg.len=1,cex=0.8)# 

################################################################################
###C-index
#options(download.file.method = 'libcurl')
#options(url.method='libcurl')
#if(!require("survcomp")) BiocManager::install("survcomp",update = F,ask = F)
library(survcomp)
library(survival)
#############################TR-DFS##################################################3
rt=read.table("nomo.score.txt",header=T,sep="\t",check.names=F,row.names=1)  

fit1 =coxph(Surv(A1_OS,A2_Event) ~ riskScore+A4_N+A3_T+A17_Age, data = rt)
cindex1 =concordance.index(predict(fit1),surv.time =rt$A1_OS,surv.event=rt$A2_Event,method ="noether")
cindex1$c.index
cindex1$lower
cindex1$upper

################################################################################
###DCA
################DFS##############
data.set <- rt
attach(data.set)
source("stdca.R")
library(survival)
Srv = Surv(data.set$A1_OS, data.set$A2_Event)
Nomogram <- coxph(Srv ~ riskScore+A4_N+A3_T+A17_Age, data=data.set)
model1 <- coxph(Srv ~ A4_N+A3_T+A17_Age, data=data.set)

cutoff <- 5

data.set$Nomogram <- c(1 - (summary(survfit(Nomogram,newdata=data.set), 
                                    times=cutoff)$surv))
data.set$model1 <- c(1 - (summary(survfit(model1,newdata=data.set), 
                                  times=cutoff)$surv))

stdca(data=data.set, outcome="A2_Event", ttoutcome="A1_OS", 
      timepoint=cutoff,predictors=c("Nomogram","model1"), 
      xstop=0.76, smooth=T)

################################################################################
################################################################################
################################################################################
###NRI IDI
##############################TR-DFS####################################
#install.packages("nricens")
library(nricens)
library(PredictABEL)

dev=read.table("nomo.score.txt",header=T,sep="\t",check.names=F,row.names=1)  

mstd=coxph(Surv(dev$A1_OS,dev$A2_Event==1)~A4_N+A3_T+A17_Age,x=TRUE,data=dev)
mnew=coxph(Surv(dev$A1_OS,dev$A2_Event==1)~riskScore+A4_N+A3_T+A17_Age,x=TRUE,data=dev)

cutoff=5
p.std= get.risk.coxph(mstd, t0=cutoff)
p.new= get.risk.coxph(mnew, t0=cutoff)
#################other##############################
#nricens(mdl.std = mstd, mdl.new = mnew,t0 = cutoff, updown = 'diff',cut = 0.05, niter = 1000)
pstd <- p.std
pnew <- p.new
reclassification(data = dev,cOutcome = 2,
                 predrisk1 = pstd,predrisk2 = pnew,
                 cutoff = c(0,0.30,0.60,1))

###IOBR
####################################################################

rm(list = ls())
gc()

library(limma)
library(dplyr)

load("STAD_TPM_Clinical.Rdata")

eset_stad <- TPM_T

cibersort<-deconvo_tme(eset = eset_stad,method = "cibersort",arrays = FALSE, perm = 200 )
xcell<-deconvo_tme(eset = eset_stad, method = "xcell",arrays = FALSE)
estimate<-deconvo_tme(eset = eset_stad, method = "estimate")
estimate$ID = gsub("-",".",estimate$ID)

tme_combine<-cibersort %>% 
  inner_join(.,xcell,by     = "ID") %>%
  inner_join(.,estimate,by  = "ID")  

save(cibersort,xcell,estimate,tme_combine,file="Immune_Cell.Rdata")

cib <- cibersort[,1:23]
colnames(cib) <- gsub("_CIBERSORT","",colnames(cib))
colnames(cib) <- gsub("_"," ",colnames(cib))
rownames(cib) <- cib$ID

xcell <- xcell
colnames(xcell) <- gsub("_xCell","",colnames(xcell))
colnames(xcell) <- gsub("_"," ",colnames(xcell))
rownames(xcell) <- xcell$ID

est <- estimate
colnames(est) <- gsub("_estimate","",colnames(est))
colnames(est) <- gsub("_"," ",colnames(est))
rownames(est) <- est$ID

data=as.matrix(est)[,-1]
#cib quan xcell est

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)
cluster <- cluster[,8:9]

cluster$cluster[cluster$riskScore < 0.777] = "low"
cluster$cluster[cluster$riskScore > 1.288] = "high"
cluster <- cluster %>% filter(cluster == "high" | cluster == "low")
cluster <- cluster[,2:3]

sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$cluster),]
data <- data[,-c(ncol(data)-1)]
gaps=c(1, as.vector(cumsum(table(data$cluster))))
xlabels=levels(factor(data$cluster))

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
data1=t(as.matrix(data[,-ncol(data)]))
pdf("barplot.pdf",height=10,width=18)
col=rainbow(nrow(data1),s=0.7,v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data1,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
for(i in 1:length(gaps)){
  j=i+1
  rect(xleft=a1[gaps[i]], ybottom = -0.01, xright = a1[gaps[j]], ytop= -0.06, col=bioCol[i])
  text((a1[gaps[i]]+a1[gaps[j]])/2,-0.035,xlabels[i],cex=2)
}
ytick2 = cumsum(data1[,ncol(data1)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data1),col=col,pch=15,bty="n",cex=1.3)
dev.off()

data=melt(data,id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Expression")
data$Expression <- as.numeric(data$Expression)

group=levels(factor(data$cluster))
data$cluster=factor(data$cluster, levels=group)
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]

ggboxplot(data, x="Immune", y="Expression", fill="cluster",
          xlab="",
          ylab="Fraction",
          legend.title="Cluster",
          width=0.8,
          palette=bioCol) + rotate_x_text(50) +
  stat_compare_means(aes(group=cluster),
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")), 
                     label="p.signif")

###estimate
ggboxplot(data, y="Expression", fill="cluster",
          xlab="",
          ylab="Fraction",
          legend.title="Cluster",
          width=0.8,
          palette=bioCol) + facet_wrap(~Immune,ncol=4,scales="free_y") +
  stat_compare_means(aes(group=cluster),
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")), 
                     label="p.signif")

########################checkpoint############################################

library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

rm(list = ls())
gc()

load("STAD_TPM_Clinical.Rdata")

checkpoint <- data.frame(t(TPM_T)[,c("CD47","SIRPA",
                                     "HLA-A","HLA-B","HLA-C","LILRB1",
                                     "CD24","SIGLEC10")])
data <- checkpoint

cluster <- read.delim("fatty.score.txt",header = T,row.names = 1)
cluster <- cluster[,8:9]

cluster$cluster[cluster$riskScore < 0.777] = "low"
cluster$cluster[cluster$riskScore > 1.288] = "high"
cluster <- cluster %>% filter(cluster == "high" | cluster == "low")
cluster <- cluster[,2:3]

sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$cluster),]
data <- data[,-c(ncol(data)-1)]
gaps=c(1, as.vector(cumsum(table(data$cluster))))
xlabels=levels(factor(data$cluster))

data=melt(data,id.vars=c("cluster"))
colnames(data)=c("cluster", "checkpoint", "Expression")
data$Expression <- as.numeric(data$Expression)

group=levels(factor(data$cluster))
data$cluster=factor(data$cluster, levels=group)
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]

ggboxplot(data, x="checkpoint", y="Expression", fill="cluster",
          xlab="",
          ylab="Fraction",
          legend.title="Cluster",
          width=0.8,
          palette=bioCol) + rotate_x_text(50) +
  stat_compare_means(aes(group=cluster),
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1),symbols=c("***", "**", "*", "ns")), 
                     label="p.signif")

###############################################################################

rt <- read.table("macrophage_list.csv",header = T,sep = ",")
load("STAD_TPM_Clinical.Rdata")
test <- TPM_T[1:10,1:10]

data <- data.frame(t(TPM_T))
x <- data[,rt$Riskscore[1:6]]

list <- read.table("list.txt",header = F,sep = "\t")
test <- TPM_T[list$V1,]
test <- na.omit(test)
test <- data.frame(t(test))

y <- test[,1:6]
y <- test[,7:20]
y <- test[,21:34]
y <- test[,35:ncol(test)]

cor <- corr.test(y,x,method = "pearson",adjust = "BH",ci = F)
cmt<-cor$r
pmt<-cor$p

ggcorrplot(cmt,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=3,
           p.mat = pmt, insig= "blank", pch.col = "red", pch.cex = 3, tl.cex = 12)

#write.csv(pmt,file = "corall.csv")

################################################################
###CCLE

library(ggcorrplot)
library(psych)

data <- read.table("CCLE_RNAseq.txt",header = T,row.names = 1,sep = "\t")

x <- data[,c("ADH4","AKR1B1","CYP4A11","NEU2","SMPD3","ST6GALNAC3")]

list <- read.table("list.txt",header = F,sep = "\t")
data <- data.frame(t(data[,6:ncol(data)]))
test <- data[list$V1,]
test <- na.omit(test)
test <- data.frame(t(test))

y <- test[,1:6]
y <- test[,7:20]
y <- test[,21:34]
y <- test[,35:ncol(test)]

cor <- corr.test(y,x,method = "pearson",adjust = "BH",ci = F)
cmt<-cor$r
pmt<-cor$p

ggcorrplot(cmt,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=3,
           p.mat = pmt, insig= "blank", pch.col = "red", pch.cex = 3, tl.cex = 12)

################################################################################

FARP1_Expression_Public <- read.csv("ST6GALNAC3 Expression Public 22Q4.csv")

library(ggplot2)
library(ggpubr)

ggplot(FARP1_Expression_Public, 
       aes(x = reorder(`Primary.Disease`,`Expression.Public.22Q4`, FUN = median),  
           y =`Expression.Public.22Q4`,color=`Primary.Disease`)) + 
  geom_boxplot()+ 
  geom_point() + 
  theme_classic(base_size = 12)+ 
  rotate_x_text(45)+ 
  theme(plot.margin=unit(rep(3,4),'cm'),legend.position="none")+ 
  xlab(NULL)+ylab("FARP1 expression \nLog2(TPM+1)")+ 
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", width = 0.5,position = position_dodge(0.9))+ 
  geom_hline(yintercept = mean(FARP1_Expression_Public$Expression.Public.22Q4), lty = 2) 

data<-FARP1_Expression_Public[FARP1_Expression_Public$`Primary.Disease` == 'Stomach Adenocarcinoma',]
ggplot(data, aes(x=reorder(`Cell.Line.Name`,-`Expression.Public.22Q4`), y=`Expression.Public.22Q4`)) + 
  geom_point(aes(size=`Expression.Public.22Q4`,color=`Expression.Public.22Q4`),stat='identity') +scale_color_continuous(low='blue' ,high='red') +
  geom_segment(aes(y = mean(data$`Expression.Public.22Q4`), 
                   x = `Cell.Line.Name`, 
                   yend = `Expression.Public.22Q4`, 
                   xend = `Cell.Line.Name`), 
               color = "black") +
  theme_classic(base_size = 14) + 
  coord_flip() + 
  xlab(NULL)+ylab("ST6GALNAC3 expression")+
  geom_hline(yintercept = mean(data$`Expression.Public.22Q4`), lty = 2)

################################################################################
rt <- data.frame(t(log2TPM))[,c("ADH4", "AKR1B1", "CYP4A11", "NEU2", "SMPD3", "ST6GALNAC3")]
clin <- clinical[,c("A1_OS","A2_Event")]
data <- merge(rt,clin,by = "row.names")

data$ADH4 <- ifelse(data$ADH4 > median(data$ADH4),"High","Low")
data$AKR1B1 <- ifelse(data$AKR1B1 > median(data$AKR1B1),"High","Low")
data$CYP4A11 <- ifelse(data$CYP4A11 > median(data$CYP4A11),"High","Low")
data$NEU2 <- ifelse(data$NEU2 > median(data$NEU2),"High","Low")
data$SMPD3 <- ifelse(data$SMPD3 > median(data$SMPD3),"High","Low")
data$ST6GALNAC3 <- ifelse(data$ST6GALNAC3 > median(data$ST6GALNAC3),"High","Low")

library(ggstatsplot)
library(ggplot2)
library(dplyr)

diamonds2 <- data[,c(2:7,9)]

p <- ggbarstats(diamonds2, A2_Event, AKR1B1, palette = 'Set2')
p$data
p
