
# data
# ~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/RNAseq_raw_86/enfermos_raw.txt  
# ~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/RNAseq_raw_86/sanos_raw.txt
# ~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/miRNAseq_raw_86/mirnas_sanos_raw.txt
# ~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/miRNAseq_raw_86/mirnas_enfermos_raw.txt

# RNASeq raw data
rna_enfermos <- read.table(file = "~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/RNAseq_raw_86/enfermos_raw.txt")
rna_sanos <- read.table(file = "~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/RNAseq_raw_86/sanos_raw.txt")
identical(rownames(rna_enfermos), rownames(rna_sanos))
# if TRUE:
rna <- cbind(rna_enfermos, rna_sanos)
# round RNASeq estimate counts
rna <- round(rna, 0)
# define experiment
countData_rna <- rna
colData_rna <- matrix(c(rep('enfermos', 86), rep('sanos', 86)))
rownames(colData_rna) <- colnames(rna)
colnames(colData_rna) <- "condition"
colData_rna <- as.data.frame(colData_rna)

library("DESeq2")
# get dataset
dds_rna <-DESeqDataSetFromMatrix(countData = countData_rna, colData = colData_rna, design = ~ condition)
# prefilter
dds_rna <- dds_rna[ rowSums(counts(dds_rna)) > 1, ]
# levels
dds_rna$condition <- factor(dds_rna$condition, levels=c("enfermos","sanos"))
dds_rna$condition <- relevel(dds_rna$condition, ref= "sanos")
# differential expression 
dds_rna <- DESeq(dds_rna)
res_rna <- results(dds_rna)
res_rna
# results
resOrdered_rna <- res_rna[order(res_rna$padj),]
summary(res_rna)
sum(res_rna$padj < 0.01, na.rm=TRUE)
res01_rna <- results(dds_rna, alpha=0.01)
summary(res01_rna)
pdf("plotMA_DEseq2_rna.pdf",width=14,height=10)
plotMA(res_rna, main="DESeq2", ylim=c(-2,2))
dev.off()
resMLE_rna <- results(dds_rna, addMLE=TRUE)
head(resMLE_rna, 4)
pdf("plotMAunsh_DEseq2_rna.pdf",width=14,height=10)
plotMA(resMLE_rna, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))
dev.off()
pdf("plotcounts_rna.pdf",width=14,height=10)
plotCounts(dds_rna, gene=which.min(res_rna$padj), intgroup="condition")
dev.off()

# save results
write.table(res_rna, file="DE_seq2_results_rna.txt", sep="\t", quote=F)
write.table(res01_rna, file="DE_seq2_results_rna_adjp_01.txt", sep="\t", quote=F)

# miRNASeq raw data 
mirna_enfermos <- read.table(file = "~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/miRNAseq_raw_86/mirnas_enfermos_raw.txt")
mirna_sanos <- read.table(file = "~/TCGA_miRNA_BC/miRNAs_resultados/subredes/diffexp/miRNAseq_raw_86/mirnas_sanos_raw.txt")
identical(rownames(mirna_enfermos), rownames(mirna_sanos))
# if TRUE:
mirna <- cbind(mirna_enfermos, mirna_sanos)
# define experiment
countData_mirna <- mirna
colData_mirna <- matrix(c(rep('enfermos', 86), rep('sanos', 86)))
rownames(colData_mirna) <- colnames(mirna)
colnames(colData_mirna) <- "condition"
colData_mirna <- as.data.frame(colData_mirna)

library("DESeq2")
# get dataset
dds_mirna <-DESeqDataSetFromMatrix(countData = countData_mirna, colData = colData_mirna, design = ~ condition)
# prefilter
dds_mirna <- dds_mirna[ rowSums(counts(dds_mirna)) > 1, ]
# levels
dds_mirna$condition <- factor(dds_mirna$condition, levels=c("enfermos","sanos"))
dds_mirna$condition <- relevel(dds_mirna$condition, ref= "sanos")
# differential expression 
dds_mirna <- DESeq(dds_mirna)
res_mirna <- results(dds_mirna)
res_mirna
# results
resOrdered_mirna <- res_mirna[order(res_mirna$padj),]
summary(res_mirna)
sum(res_mirna$padj < 0.01, na.rm=TRUE)
res01_mirna <- results(dds_mirna, alpha=0.01)
summary(res01_mirna)
pdf("plotMA_DEseq2_mirna.pdf",width=14,height=10)
plotMA(res_mirna, main="DESeq2", ylim=c(-2,2))
dev.off()
resMLE_mirna <- results(dds_mirna, addMLE=TRUE)
head(resMLE_mirna, 4)
pdf("plotMAunsh_DEseq2_mirna.pdf",width=14,height=10)
plotMA(resMLE_mirna, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))
dev.off()
pdf("plotcounts_mirna.pdf",width=14,height=10)
plotCounts(dds_mirna, gene=which.min(res_mirna$padj), intgroup="condition")
dev.off()

# save results
write.table(res_mirna, file="DE_seq2_results_mirna.txt", sep="\t", quote=F)
write.table(res01_mirna, file="DE_seq2_results_mirna_adjp_01.txt", sep="\t", quote=F)


