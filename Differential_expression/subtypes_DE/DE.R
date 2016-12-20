# differential expression subtypes:
#/mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/subtipos_exp_diferencial
# R
#----------------------------------------------------------------
#----------------------------------------------------------------
# sanos:
#----------------------------------------------------------------
# genes: 	GeneID
#			sanos_raw_counts.txt
# paste GeneID sanos_raw_counts.txt  > sanos_raw_full.txt
# grep ^[^?] sanos_raw_full.txt > sanos_raw_filtered.txt
#----------------------------------------------------------------
# mirnas: 	miRNAs_sanos_maduros_filtro.txt
#----------------------------------------------------------------
#----------------------------------------------------------------
# subtipos:
#----------------------------------------------------------------
#			Subtipos_pam50robust.scale_752_RNAseq.txt
#----------------------------------------------------------------
# genes: 	raw_752.txt
# grep ^[^?] raw_752.txt > raw_rnaseq_752_filtered.txt
#----------------------------------------------------------------
# mirnas: 	miRNAs_752_maduros_filtro.txt
#----------------------------------------------------------------
#----------------------------------------------------------------

#################################################################
# Subtypes
#################################################################

subtipos <- read.table("Subtipos_pam50robust.scale_752_RNAseq.txt", header=FALSE)

subtipos_ord <- subtipos[order(subtipos[,2]),]

# basal_p <- subtipos[(subtipos[,2] == 'Basal'),]
# her2_p <- subtipos[(subtipos[,2] == 'Her2'),]
# luma_p <- subtipos[(subtipos[,2] == 'LumA'),]
# lumb_p <- subtipos[(subtipos[,2] == 'LumB'),]

#################################################################
# RNA-Seq
################################################################

# ---------------------------------------------------------------
# tumours
# ---------------------------------------------------------------

RNAseq <- read.table('raw_rnaseq_752_filtered.txt', header=T, row.names=1)

RNAseq_ord <- RNAseq[, as.vector(subtipos_ord[,1])]

# /home/diana/R/x86_64-pc-linux-gnu-library/3.3
# basal <- RNAseq[,(basal_p[,1])]
# luma <- RNAseq[,(luma_p[,1])]
# lumb <- RNAseq[,(lumb_p[,1])]
# her2 <- RNAseq[,(her2_p[,1])]

# ---------------------------------------------------------------
# controls
# ---------------------------------------------------------------

sanos <- read.table('sanos_raw_filtered.txt', header=T, row.names=1)

# ---------------------------------------------------------------
# create count and metadata data frames 'countData' & 'colData'
# ---------------------------------------------------------------

identical(rownames(RNAseq_ord), rownames(sanos))

countData_rna <- round((cbind(RNAseq_ord, sanos)), 0)

colData_rna <- matrix(c(as.vector(subtipos_ord[,2]), rep('sanos',  length(colnames(sanos)))))

rownames(colData_rna) <- colnames(countData_rna)

colnames(colData_rna) <- "condition"

colData_rna <- as.data.frame(colData_rna)

# ---------------------------------------------------------------
# DIferential expression RNASeq
# ---------------------------------------------------------------

library("DESeq2")

dds_rna <-DESeqDataSetFromMatrix(countData = countData_rna, colData = colData_rna, design = ~ condition)

dds_rna <- dds_rna[ rowSums(counts(dds_rna)) > 1, ]

dds_rna$condition <- relevel(dds_rna$condition, ref= "sanos")

dds_rna <- DESeq(dds_rna)
res_rna <- results(dds_rna)
res_rna

# ---------------------------------------------------------------
# Differential expression RNASeq Subtypes
# ---------------------------------------------------------------

res_basal <- results(dds_rna, contrast=c("condition","Basal","sanos"))
write.table(res_basal, file="DE_basal.txt", sep="\t", quote=F, col.names=NA, row.names=T)

res_her2 <- results(dds_rna, contrast=c("condition","Her2","sanos"))
write.table(res_her2, file="DE_her2.txt", sep="\t", quote=F, col.names=NA, row.names=T)

res_luma <- results(dds_rna, contrast=c("condition","LumA","sanos"))
write.table(res_luma, file="DE_luma.txt", sep="\t", quote=F, col.names=NA, row.names=T)

res_lumb <- results(dds_rna, contrast=c("condition","LumB","sanos"))
write.table(res_lumb, file="DE_lumb.txt", sep="\t", quote=F, col.names=NA, row.names=T)


#################################################################
# micro-RNA-Seq
################################################################


# ---------------------------------------------------------------
# tumours
# ---------------------------------------------------------------

miRNAseq <- read.table('miRNAs_752_maduros_filtro.txt', header=T, row.names=1)

miRNAseq_ord <- miRNAseq[, as.vector(subtipos_ord[,1])]

# ---------------------------------------------------------------
# controls
# ---------------------------------------------------------------

sanos_mir <- read.table('miRNAs_sanos_maduros_filtro.txt', header=T, row.names=1)

# ---------------------------------------------------------------
# create count and metadata data frames 'countData' & 'colData'
# ---------------------------------------------------------------



miRNAseq_ord <- miRNAseq_ord[(intersect(rownames(sanos_mir), rownames(miRNAseq_ord))),]

identical(colnames(miRNAseq_ord), as.vector(subtipos_ord[,1]))

identical(rownames(miRNAseq_ord), rownames(sanos_mir))

countData_mirna <- round((cbind(miRNAseq_ord, sanos_mir)), 0)

colData_mirna <- matrix(c(as.vector(subtipos_ord[,2]), rep('sanos',  length(colnames(sanos_mir)))))

rownames(colData_mirna) <- colnames(countData_mirna)

colnames(colData_mirna) <- "condition"

colData_mirna <- as.data.frame(colData_mirna)

# ---------------------------------------------------------------
# Diferential expression miRNASeq
# ---------------------------------------------------------------

library("DESeq2")

dds_mirna <-DESeqDataSetFromMatrix(countData = countData_mirna, colData = colData_mirna, design = ~ condition)

dds_mirna <- dds_mirna[ rowSums(counts(dds_mirna)) > 1, ]

dds_mirna$condition <- relevel(dds_mirna$condition, ref= "sanos")

dds_mirna <- DESeq(dds_mirna)
res_mirna <- results(dds_mirna)
res_mirna

# ---------------------------------------------------------------
# Differential expression miRNASeq Subtypes
# ---------------------------------------------------------------

res_mir_basal <- results(dds_mirna, contrast=c("condition","Basal","sanos"))
write.table(res_mir_basal, file="DE_mir_basal.txt", sep="\t", quote=F, col.names=NA, row.names=T)

res_mir_her2 <- results(dds_mirna, contrast=c("condition","Her2","sanos"))
write.table(res_mir_her2, file="DE_mir_her2.txt", sep="\t", quote=F, col.names=NA, row.names=T)

res_mir_luma <- results(dds_mirna, contrast=c("condition","LumA","sanos"))
write.table(res_mir_luma, file="DE_mir_luma.txt", sep="\t", quote=F, col.names=NA, row.names=T)

res_mir_lumb <- results(dds_mirna, contrast=c("condition","LumB","sanos"))
write.table(res_mir_lumb, file="DE_mir_lumb.txt", sep="\t", quote=F, col.names=NA, row.names=T)

save.image("mir_DE_subtypes.RData")
