

# read raw counts matrix and Gene IDs
RNAseq <- read.table(file = "raw_counts.txt", header=T, sep = '\t')
RNAseq <- RNAseq[30:20531, ]
ID <- scan(file="GeneID", skip=30, what="")
row.names(RNAseq) <- ID
#RNAseq <- as.matrix(RNAseq)

# obtain gene list of genes used for the mi network analysis
setwd("/home/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/")
enfermos <- read.table(file = "enfermos.86.txt", header = T, sep = '\t')
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')
enfermos <- as.vector(enfermos[ ,1])
sanos <- as.vector(sanos[ ,1])
genes <- sort(intersect(enfermos, sanos))

# filter raw count matrix by gene list
RNAseq <- RNAseq[genes, ]

 round raw count matrix values
RNAseq <- round(RNAseq, 0)

# define count table and metadata
countData <- RNAseq
colData <- matrix(c(rep('enfermos', 86), rep('sanos', 86)))
rownames(colData) <- colnames(RNAseq)
colnames(colData) <- "condition"

colData <- as.data.frame(colData)
#countData <- as.data.frame(countData)

library("DESeq2")

dds<-DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds$condition <- relevel(dds$condition, ref= "sanos")

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)


pdf("plotMA_DEseq2_RNAseq.pdf",width=14,height=10)
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)

pdf("plotMAunsh_DEseq2_RNAseq.pdf",width=14,height=10)
plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))
dev.off()

pdf("plotcounts_RNAseq.pdf",width=14,height=10)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

# save results
write.table(res, file="DE_seq2_results.txt", sep="\t", quote=F)