# EBseq
# the matrix must contain RSEM raw counts: Genes x Samples (all groups in the same matrix)
# paste enfermos_raw_counts.txt sanos_raw_counts.txt > raw_counts_RNAseq.txt 

library(EBSeq)

# read matrix and Gene IDs
RNAseq <- read.table(file = "raw_counts.txt", header=T, sep = '\t')
RNAseq <- RNAseq[30:20531, ]
ID <- scan(file="GeneID", skip=30, what="")
row.names(RNAseq) <- ID
RNAseq <- as.matrix(RNAseq)

setwd("/home/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/")

# read expression matrices that were used for mi computation
enfermos <- read.table(file = "enfermos.86.txt", header = T, sep = '\t')
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')

# create a vector with gene IDs from both matrices
enfermos <- as.vector(enfermos[ ,1])
sanos <- as.vector(sanos[ ,1])


# gene list from genes in common
genes <- sort(intersect(enfermos, sanos))

# filter matrix by gene list 'genes'
RNAseq <- RNAseq[genes, ]

# define factors, library size factor: using Upper-Quartile normalization as TCGA (as DEseq: Sizes=MedianNorm(RNAseq))
myfactors <- as.factor(c( rep('enfermos', 86), rep('sanos', 86)))
str(RNAseq)
Sizes=QuantileNorm(RNAseq,.75)

# run EB Test
EBOut=EBTest(Data=RNAseq, Conditions=myfactors, sizeFactors=Sizes, maxround=7)

# check for convergence < 0.01
EBOut$Alpha
EBOut$Beta
EBOut$P

# obtain results
EBDERes=GetDEResults(EBOut, FDR=0.05)
str(EBDERes$DEfound)

# save diferential expressed genes 
write.table(x=EBDERes$DEfound, file="EBseq_dif_exp_FDR0.05_filtrogenes", sep="\t", quote=F)

# run Fold Change between two conditions
GeneFC=PostFC(EBOut)
str(GeneFC)

pdf("post_vs_raw_fc_filtrogenes.pdf",width=14,height=10)
PlotPostVsRawFC(EBOut,GeneFC)
dev.off()

# volcano plot
## from Michael S. Chimenti, PhD http://www.michaelchimenti.com/2015/08/create-a-volcano-plot-on-ebseq-output/
volc_df = data.frame(EBOut$PPDE, GeneFC$PostFC)
pdf("volcano_ebseq_rnaseq_filtrogenes.pdf",width=14,height=10)
with(volc_df, plot(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, main="Volcano Plot EBSeq", xlim=c(-10,10)))
abline(h=0.95)
with(subset(volc_df, EBOut.PPDE > 0.95 & abs(log2(GeneFC.PostFC)) < 1.5), points(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, col="orange")) 
with(subset(volc_df, EBOut.PPDE > 0.95 & abs(log2(GeneFC.PostFC)) > 1.5), points(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, col="red"))
dev.off()

# diferential expression
volc_df[, 2] <- log2(volc_df[2])
dif_exp <- volc_df[EBDERes$DEfound, ]
write.table(dif_exp, quote = F, sep = "\t", row.names = T, col.names = NA, file = "EBseq_dif_exp_only_filtrogenes.txt")
write.table(volc_df, quote = F, sep = "\t", row.names = T, col.names = NA, file = "EBseq_PPDE_PFC_filtrogenes.txt")

matrix_norm <- GetNormalizedMat(Data=RNAseq, Sizes=Sizes)
write.table(matrix_norm, file="Matriz_normalizada_RNAseq_filtrogenes.txt", sep="\t", quote=F)

#################################################RESULTS#########################################################################
# EBOut=EBTest(Data=RNAseq, Conditions=myfactors, sizeFactors=Sizes, maxround
=8)                                                                             
# iteration 1 done 
# time 287.19 
# iteration 2 done 
# time 117.4 
# iteration 3 done 
# time 110.98 
# iteration 4 done 
# time 110.64 
# iteration 5 done 
# time 113.96 
# iteration 6 done 
# time 100.53 
# iteration 7 done 
# time 100.78 
# iteration 8 done 
# time 100.72 
# EBOut1$Alpha
#           [,1]
# iter1 1.300362
# iter2 1.310479
# iter3 1.323406
# iter4 1.333003
# iter5 1.352608
# iter6 1.352608
# iter7 1.352608
# iter8 1.352608
# EBOut1$Beta
#             Ng1
# iter1  97.43924
# iter2 101.06623
# iter3 101.84036
# iter4 103.26316
# iter5 104.79683
# iter6 105.45181
# iter7 105.45181
# iter8 105.45181
# EBOut1$P
#            [,1]
# iter1 0.6044040
# iter2 0.7295917
# iter3 0.7654530
# iter4 0.7796165
# iter5 0.7841911
# iter6 0.7841911
# iter7 0.7841911
# iter8 0.7841911
# with MedianNorm
# length(EBDERes1$DEfound)
# [1] 3322
# with upper quantile
# length(EBDERes1$DEfound)
# [1] 3253


########################################## analysis with all genes ###########################################################

# library(EBSeq)
# read matrix and Gene IDs
# RNAseq <- read.table(file = "raw_counts.txt", header=T, sep = '\t')
# RNAseq <- RNAseq[30:20531, ]
# ID <- scan(file="GeneID", skip=30, what="")
# row.names(RNAseq) <- ID
# RNAseq <- as.matrix(RNAseq)
# define factors, library size factor: using Upper-Quartile normalization as TCGA (as DEseq: Sizes=MedianNorm(RNAseq))
# myfactors <- as.factor(c( rep('enfermos', 86), rep('sanos', 86)))
# str(RNAseq)
# Sizes=QuantileNorm(RNAseq,.75)
# run EB Test
# EBOut=EBTest(Data=RNAseq, Conditions=myfactors, sizeFactors=Sizes, maxround=5)
# EBOut$Alpha
# EBOut$Beta
# EBOut$P
# EBDERes=GetDEResults(EBOut, FDR=0.05)
# str(EBDERes$DEfound)
# length(EBDERes$DEfound)
# save diferential expressed genes 
# write.table(x=EBDERes$DEfound, file="EBseq_dif_exp_FDR0.05", sep="\t", quote=F)
# run Fold Change between two conditions
# GeneFC=PostFC(EBOut)
# str(GeneFC)
# pdf("post_vs_raw_fc.pdf",width=14,height=10)
# PlotPostVsRawFC(EBOut,GeneFC)
# dev.off()
# volcano plot
## from Michael S. Chimenti, PhD http://www.michaelchimenti.com/2015/08/create-a-volcano-plot-on-ebseq-output/
# volc_df = data.frame(EBOut$PPDE, GeneFC$PostFC)
# pdf("volcano_ebseq_rnaseq.pdf",width=14,height=10)
# with(volc_df, plot(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, main="Volcano Plot EBSeq", xlim=c(-12,12)))
# abline(h=0.95)
# with(subset(volc_df, EBOut.PPDE > 0.95 & abs(log2(GeneFC.PostFC)) < 1.5), points(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, col="orange")) 
# with(subset(volc_df, EBOut.PPDE > 0.95 & abs(log2(GeneFC.PostFC)) > 1.5), points(log2(GeneFC.PostFC), EBOut.PPDE, pch=20, col="red"))
# dev.off()
# diferential expression
# volc_df[, 2] <- log2(volc_df[2])
# dif_exp <- volc_df[EBDERes$DEfound, ]
# write.table(dif_exp, quote = F, sep = "\t", row.names = T, col.names = NA, file = "EBseq_dif_exp_only.txt")
# write.table(volc_df, quote = F, sep = "\t", row.names = T, col.names = NA, file = "EBseq_PPDE_PFC.txt")
# matrix_norm <- GetNormalizedMat(Data=RNAseq, Sizes=Sizes)
# write.table(matrix_norm, file="Matriz_normalizada_RNAseq", sep="\t", quote=F)
#################################################RESULTS#########################################################################
# EBOut=EBTest(Data=RNAseq, Conditions=myfactors, sizeFactors=Sizes, maxround=5) 
# Removing transcripts with 100 th quantile < = 0 
# 20026 transcripts will be tested
# iteration 1 done 
# time 175.86 
# iteration 2 done 
# time 95.11 
# iteration 3 done 
# time 90.44 
# iteration 4 done 
# time 76.75
# iteration 5 done 
# time 76.72 
#  EBOut$Alpha
#            [,1]
# iter1 0.3124202
# iter2 0.3147973
# iter3 0.3129595
# iter4 0.3129595
# iter5 0.3129595
#  EBOut$Beta
#            Ng1
# iter1 1.368995
# iter2 1.475926
# iter3 1.474103
# iter4 1.474103
# iter5 1.474103
#  EBOut$P
#            [,1]
# iter1 0.6508818
# iter2 0.6916856
# iter3 0.7012750
# iter4 0.7012750
# iter5 0.7012750
# with MedianNorm
# length(EBDERes$DEfound)
# [1] 4722
# with upper quantile
# length(EBDERes$DEfound)
# [1] 4682 

###############################################################
