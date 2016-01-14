

# read raw counts matrix and Gene IDs
RNAseq1 <- read.table(file = "enfermos_raw_752.txt", header=T, sep = '\t')
RNAseq2 <- read.table(file = "sanos_raw_counts.txt", header=T, sep = '\t')
RNAseq1 <- RNAseq1[30:20531, ]
RNAseq2 <- RNAseq2[30:20531, ]

RNAseq <- cbind(RNAseq1, RNAseq2)

ID <- scan(file="GeneID", skip=30, what="")
row.names(RNAseq) <- ID

# obtain gene list of genes used for the mi network analysis
enfermos <- read.table("RNAseq_752_norm.rsem.aracne.txt", header=T, sep="\t")
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')
enfermos <- as.vector(enfermos[ ,1])
sanos <- as.vector(sanos[ ,1])
genes <- sort(intersect(enfermos, sanos))

# length(genes)
#[1] 15219

# filter raw count matrix by gene list
RNAseq <- RNAseq[genes, ]

# round raw count matrix values
RNAseq <- round(RNAseq, 0)

# define count table and metadata
countData <- RNAseq
subtype <- read.table("Subtipos_pam50.scale_752_RNAseq.txt", sep="\t", row.names=1)
colnames(subtype) <- "condition"
controles1 <- colnames(RNAseq2)
controles2 <- rep("Control", 86)
controles <- as.data.frame(controles2, row.names=controles1)
colnames(controles) <- "condition"
colData <- rbind(subtype, controles)
col1 <- matrix(c(rep('enfermos', 752), rep('sanos', 86)))
colnames(col1) <- "condition2"
colData <- cbind(colData, col1)

library("DESeq2")

dds<-DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition2)
dds$condition2 <- relevel(dds$condition2, ref= "sanos")

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)


pdf("plotMA_DEseq2_RNAseq752.pdf",width=14,height=10)
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)

pdf("plotMAunsh_DEseq2_RNAseq752.pdf",width=14,height=10)
plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))
dev.off()

pdf("plotcounts_RNAseq752.pdf",width=14,height=10)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

# save results
write.table(res, file="DE_seq2_results752.txt", sep="\t", quote=F)

##############################################
library("BiocParallel")
register(MulticoreParam(10))

### subtypes

dds1 <-DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1$condition <- relevel(dds1$condition, ref= "Control")

dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1

res_luma <- results(dds1, contrast=c("condition","LumA","Control"))
write.table(res_luma, file="DE_luma_control.txt", sep="\t", quote=F)
res_lumavslumb <-results(dds1, contrast=c("condition","LumA","LumB"))
write.table(res_lumavslumb, file="DE_luma_lumb.txt", sep="\t", quote=F)
res_lumavsher2 <-results(dds1, contrast=c("condition","LumA","Her2"))
write.table(res_lumavsher2, file="DE_luma_her2.txt", sep="\t", quote=F)
res_lumavsbasal <-results(dds1, contrast=c("condition","LumA","Basal"))
write.table(res_lumavsbasal, file="DE_luma_basal.txt", sep="\t", quote=F)
res_lumavsnormal <-results(dds1, contrast=c("condition","LumA","Normal"))
write.table(res_lumavsnormal, file="DE_luma_normal.txt", sep="\t", quote=F)

res_lumb <-results(dds1, contrast=c("condition","LumB","Control"))
write.table(res_lumb, file="DE_lumb_control.txt", sep="\t", quote=F)
res_lumbvsluma <-results(dds1, contrast=c("condition","LumB","LumA"))
write.table(res_lumbvsluma, file="DE_lumb_luma.txt", sep="\t", quote=F)
res_lumbvsher2 <-results(dds1, contrast=c("condition","LumB","Her2"))
write.table(res_lumbvsher2, file="DE_lumb_her2.txt", sep="\t", quote=F)
res_lumbvsbasal <-results(dds1, contrast=c("condition","LumB","Basal"))
write.table(res_lumbvsbasal, file="DE_lumb_basal.txt", sep="\t", quote=F)
res_lumbvsnormal <-results(dds1, contrast=c("condition","LumB","Normal"))
write.table(res_lumbvsnormal, file="DE_lumb_normal.txt", sep="\t", quote=F)

res_her2 <-results(dds1, contrast=c("condition","Her2","Control"))
write.table(res_her2, file="DE_her2_control.txt", sep="\t", quote=F)
res_her2vsluma <-results(dds1, contrast=c("condition","Her2","LumA"))
write.table(res_her2vsluma, file="DE_her2_luma.txt", sep="\t", quote=F)
res_her2vslumb <-results(dds1, contrast=c("condition","Her2","LumB"))
write.table(res_her2vslumb, file="DE_her2_lumb.txt", sep="\t", quote=F)
res_her2vsbasal <-results(dds1, contrast=c("condition","Her2","Basal"))
write.table(res_her2vsbasal, file="DE_her2_basal.txt", sep="\t", quote=F)
res_her2vsnormal <-results(dds1, contrast=c("condition","Her2","Normal"))
write.table(res_her2vsnormal, file="DE_her2_normal.txt", sep="\t", quote=F)

res_basal <-results(dds1, contrast=c("condition","Basal","Control"))
write.table(res_basal, file="DE_basal_control.txt", sep="\t", quote=F)
res_basalvsluma <-results(dds1, contrast=c("condition","Basal","LumA"))
write.table(res_basalvsluma, file="DE_basal_luma.txt", sep="\t", quote=F)
res_basalvslumb <-results(dds1, contrast=c("condition","Basal","LumB"))
write.table(res_basalvslumb, file="DE_basal_lumb.txt", sep="\t", quote=F)
res_basalvsher2 <-results(dds1, contrast=c("condition","Basal","Her2"))
write.table(res_basalvsher2, file="DE_basal_her2.txt", sep="\t", quote=F)
res_basalvsnormal <-results(dds1, contrast=c("condition","Basal","Normal"))
write.table(res_basalvsnormal, file="DE_basal_normal.txt", sep="\t", quote=F)

res_normal <-results(dds1, contrast=c("condition","Normal","Control"))
write.table(res_normal, file="DE_normal_control.txt", sep="\t", quote=F)
res_normalvsluma <-results(dds1, contrast=c("condition","Normal","LumA"))
write.table(res_normalvsluma, file="DE_normal_luma.txt", sep="\t", quote=F)
res_normalvslumb <-results(dds1, contrast=c("condition","Normal","LumB"))
write.table(res_normalvslumb, file="DE_normal_lumb.txt", sep="\t", quote=F)
res_normalvsher2 <-results(dds1, contrast=c("condition","Normal","Her2"))
write.table(res_normalvsher2, file="DE_normal_her2.txt", sep="\t", quote=F)
res_normalvsbasal <-results(dds1, contrast=c("condition","Normal","Basal"))
write.table(res_normalvsbasal, file="DE_normal_basal.txt", sep="\t", quote=F)



pdf("plotMA_DEseq2_RNAseq752.pdf",width=14,height=10)
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)

pdf("plotMAunsh_DEseq2_RNAseq752.pdf",width=14,height=10)
plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))
dev.off()

pdf("plotcounts_RNAseq752.pdf",width=14,height=10)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

# save results
write.table(res, file="DE_seq2_results752.txt", sep="\t", quote=F)

