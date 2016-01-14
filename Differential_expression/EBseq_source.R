
library(EBSeq)
RNAseq <- read.table(file = "raw_counts.txt", header=T, sep = '\t')
RNAseq <- RNAseq[30:20531, ]
ID <- scan(file="GeneID", skip=30, what="")
row.names(RNAseq) <- ID
RNAseq <- as.matrix(RNAseq)
setwd("/home/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/")
enfermos <- read.table(file = "enfermos.86.txt", header = T, sep = '\t')
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')
enfermos <- as.vector(enfermos[ ,1])
sanos <- as.vector(sanos[ ,1])
genes <- sort(intersect(enfermos, sanos))
RNAseq <- RNAseq[genes, ]
myfactors <- as.factor(c( rep('enfermos', 86), rep('sanos', 86)))
str(RNAseq)
Sizes=QuantileNorm(RNAseq,.75)
EBOut=EBTest(Data=RNAseq, Conditions=myfactors, sizeFactors=Sizes, maxround=7)
EBOut$Alpha
EBOut$Beta
EBOut$P
EBDERes=GetDEResults(EBOut, FDR=0.05)
str(EBDERes$DEfound)
length(EBDERes$DEfound)
write.table(x=EBDERes$DEfound, file="EBseq_dif_exp_FDR0.05_filtrogenes", sep="\t", quote=F)
GeneFC=PostFC(EBOut)
str(GeneFC)
pdf("post_vs_raw_fc_filtrogenes.pdf",width=14,height=10)
PlotPostVsRawFC(EBOut,GeneFC)
dev.off()
volc_df = data.frame(EBOut$PPDE, GeneFC$PostFC)
pdf("volcano_ebseq_rnaseq_filtrogenes.pdf",width=14,height=10)
with(volc_df, plot(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, main="Volcano Plot EBSeq", xlim=c(-10,10)))
abline(h=0.95)
with(subset(volc_df, EBOut.PPDE > 0.95 & abs(log2(GeneFC.PostFC)) < 1.5), points(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, col="orange")) 
with(subset(volc_df, EBOut.PPDE > 0.95 & abs(log2(GeneFC.PostFC)) > 1.5), points(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, col="red"))
dev.off()
volc_df[, 2] <- log2(volc_df[2])
dif_exp <- volc_d[EBDERes$DEfound, ]
write.table(dif_exp, quote = F, sep = "\t", row.names = T, col.names = NA, file = "EBseq_dif_exp_only_filtrogenes.txt")
write.table(volc_df, quote = F, sep = "\t", row.names = T, col.names = NA, file = "EBseq_PPDE_PFC_filtrogenes.txt")
matrix_norm <- GetNormalizedMat(Data=RNAseq, Sizes=Sizes)
write.table(matrix_norm, file="Matriz_normalizada_RNAseq_filtrogenes.txt", sep="\t", quote=F)


# awk -F'\t' '$3 > 1.5 {print $1"\t"$2"\t"$3}' EBseq_dif_exp_only_filtrogenes.txt | less -S
# awk -F'\t' '$3 > 1.5 {print $1"\t"$2"\t"$3}' EBseq_dif_exp_only_filtrogenes.txt | wc -l
# 459
# awk -F'\t' '$3 < -1.5 {print $1"\t"$2"\t"$3}' EBseq_dif_exp_only_filtrogenes.txt | less -S
# awk -F'\t' '$3 < -1.5 {print $1"\t"$2"\t"$3}' EBseq_dif_exp_only_filtrogenes.txt | wc -l
# 444
#
# 901 dif_exp
