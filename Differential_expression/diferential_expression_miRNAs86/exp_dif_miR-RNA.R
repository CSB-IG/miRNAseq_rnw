miRNA exp dif

RNAseq <- read.table(file = "raw_counts.txt", header=T, sep = '\t')
RNAseq <- RNAseq[30:20531, ]
ID <- scan(file="GeneID", skip=30, what="")
row.names(RNAseq) <- ID

genesRNA <- scan(file="lista.DE.mir.txt", what="")
# awk '{print $1}' enfermos_RNAseq.adj.txt > lista.DE.mir.txt

# matriz 86 enfermos/sanos con 15136 genes
RNAseq <- RNAseq[genesRNA, ]
RNAseq <- round(RNAseq, 0)

# matriz mirnas

enfermos <- read.table("miRNAs_enfermos_maduros_filtro.txt", sep="\t", header=T, row.names=1)
sanos <- read.table("miRNAs_sanos_maduros_filtro.txt", sep="\t", header=T, row.names=1)

identical(rownames(enfermos), rownames(sanos))

# si TRUE:

MIRNA <- cbind(enfermos, sanos)

# MIRNA y RNAseq dim() iguales

identical(colnames(RNAseq), colnames(MIRNA))

# si TRUE:

countData <- rbind(RNAseq, MIRNA)
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


# save results
write.table(res, file="DE_seq2_results_mirnas86.txt", sep="\t", quote=F)
