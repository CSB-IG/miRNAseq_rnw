# differential expression DESeq2:
#----------------------------------------------------------------
# raw_expression_rnaseq_tumour.txt: a matrix containing raw cases RNA-Seq expression data
# raw_expression_rnaseq_control.txt: a matrix containing raw control RNA-Seq expression data
# raw_expression_mature_mirnaseq_tumour.txt: a matrix containing raw tumour miRNA-Seq expression data
# raw_expression_mature_mirnaseq_tumour.txt: a matrix containing raw control miRNA-Seq expression data
#-----------------------------------------------------------------
# Tab delimited non quoted .txt file matrices with the following format:
# Name	Sample1	Sample2	Sample3
# Gene1	Value	Value	Value
# Gene2	Value	Value	Value
# Gene3	Value	Value	Value
# ----------------------------------------------------------------

#################################################################
# RNA-Seq
################################################################

# ---------------------------------------------------------------
# tumours
# ---------------------------------------------------------------

RNAseq_tumour <- read.table('raw_expression_rnaseq_tumour.txt', header=T, row.names=1)

# ---------------------------------------------------------------
# controls
# ---------------------------------------------------------------

RNAseq_control <- read.table('raw_expression_rnaseq_control.txt', header=T, row.names=1)

# ---------------------------------------------------------------
# create count and metadata data frames 'countData' & 'colData'
# ---------------------------------------------------------------

if (identical(rownames(RNAseq_tumour), rownames(RNAseq_control)) == TRUE){
	countData_rna <- round((cbind(RNAseq_tumour, RNAseq_control)), 0)

} else {RNAseq_tumour <- RNAseq_tumour[(intersect(rownames(RNAseq_control), rownames(RNAseq_tumour))),]
	countData_rna <- round((cbind(RNAseq_tumour, RNAseq_control)), 0)
}

colData_rna <- matrix(c(rep('tumour',  length(colnames(RNAseq_tumour))), rep('control',  length(colnames(RNAseq_control)))))
rownames(colData_rna) <- colnames(countData_rna)
colnames(colData_rna) <- "condition"
colData_rna <- as.data.frame(colData_rna)

# ---------------------------------------------------------------
# DIferential expression RNASeq
# ---------------------------------------------------------------

library("DESeq2")

dds_rna <-DESeqDataSetFromMatrix(countData = countData_rna, colData = colData_rna, design = ~ condition)

dds_rna <- dds_rna[ rowSums(counts(dds_rna)) > 1, ]

dds_rna$condition <- relevel(dds_rna$condition, ref= "control")

dds_rna <- DESeq(dds_rna)

# ---------------------------------------------------------------
# Results RNA-Seq
# ---------------------------------------------------------------

res_rna <- results(dds_rna)
res_rna

# ---------------------------------------------------------------
# Save Results RNASeq
# ---------------------------------------------------------------

write.table(res_rna, file="DE_RNAseq.txt", sep="\t", quote=F, col.names=NA, row.names=T)

#################################################################
# micro-RNA-Seq
################################################################


# ---------------------------------------------------------------
# tumours
# ---------------------------------------------------------------

miRNAseq_tumour <- read.table('raw_expression_mature_mirnaseq_tumour.txt', header=T, row.names=1)

# ---------------------------------------------------------------
# controls
# ---------------------------------------------------------------

control_mir <- read.table('raw_expression_mature_mirnaseq_control.txt', header=T, row.names=1)

# ---------------------------------------------------------------
# create count and metadata data frames 'countData' & 'colData'
# ---------------------------------------------------------------



miRNAseq_tumour <- miRNAseq_tumour[(intersect(rownames(miRNAseq_control), rownames(miRNAseq_tumour))),]

identical(rownames(miRNAseq_tumour), rownames(miRNAseq_control))

countData_mirna <- round((cbind(miRNAseq_tumour, miRNAseq_control)), 0)

colData_mirna <- matrix(c(rep('tumour',  length(colnames(miRNAseq_tumour))), rep('control',  length(colnames(miRNAseq_control)))))

rownames(colData_mirna) <- colnames(countData_mirna)

colnames(colData_mirna) <- "condition"

colData_mirna <- as.data.frame(colData_mirna)

# ---------------------------------------------------------------
# Diferential expression miRNASeq
# ---------------------------------------------------------------

library("DESeq2")

dds_mirna <-DESeqDataSetFromMatrix(countData = countData_mirna, colData = colData_mirna, design = ~ condition)

dds_mirna <- dds_mirna[ rowSums(counts(dds_mirna)) > 1, ]

dds_mirna$condition <- relevel(dds_mirna$condition, ref= "control")

dds_mirna <- DESeq(dds_mirna)

# ---------------------------------------------------------------
# Differential expression miRNASeq
# ---------------------------------------------------------------

res_mirna <- results(dds_mirna)
res_mirna

# ---------------------------------------------------------------
# Save Results miRNASeq
# ---------------------------------------------------------------

write.table(res_mirna, file="DE_mir.txt", sep="\t", quote=F, col.names=NA, row.names=T)
