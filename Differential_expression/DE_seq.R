DEseq

# read raw counts matrix and Gene IDs
RNAseq <- read.table(file = "raw_counts.txt", header=T, sep = '\t')
RNAseq <- RNAseq[30:20531, ]
ID <- scan(file="GeneID", skip=30, what="")
row.names(RNAseq) <- ID
RNAseq <- as.matrix(RNAseq)

# obtain gene list of genes used for the mi network analysis
setwd("/home/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/")
enfermos <- read.table(file = "enfermos.86.txt", header = T, sep = '\t')
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')
enfermos <- as.vector(enfermos[ ,1])
sanos <- as.vector(sanos[ ,1])
genes <- sort(intersect(enfermos, sanos))

# filter raw count matrix by gene list
RNAseq <- RNAseq[genes, ]

# round raw count matrix values
RNAseq <- round(RNAseq, 0)

# define count table and metadata
countTable = RNAseq
condition = factor(c( rep('enfermos', 86), rep('sanos', 86)))

library( "DESeq" )

# create DEseq input
cds = newCountDataSet( countTable, condition )

# normalization
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )

# variance estimation
cds = estimateDispersions( cds )
str( fitInfo(cds) )
pdf("DE_seq_disp.pdf",width=14,height=10)
plotDispEsts( cds )
dev.off()
head( fData(cds) )

# differential expression anlysis: cancer(B)/control(A)
res = nbinomTest( cds, "sanos", "enfermos")
head(res)

pdf("plotMA_RNAseq.pdf",width=14,height=10)
plotMA(res)
dev.off()

pdf("histbinom_RNAseq.pdf",width=14,height=10)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

resSig = res[ res$padj < 0.1, ]

head( resSig[ order(resSig$pval), ] )

downreagulated
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )

up regulated 
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )

# save results
write.table(res, file="DE_seq_results.txt", sep="\t", quote=F, rownames=F)


##### source 

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
RNAseq <- round(RNAseq, 0)
countTable = RNAseq
condition = factor(c( rep('enfermos', 86), rep('sanos', 86)))
library( "DESeq" )
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
str( fitInfo(cds) )
pdf("DE_seq_disp.pdf",width=14,height=10)
plotDispEsts( cds )
dev.off()
head( fData(cds) )
res = nbinomTest( cds, "sanos", "enfermos")
head(res)
pdf("plotMA_RNAseq.pdf",width=14,height=10)
plotMA(res)
dev.off()
pdf("histbinom_RNAseq.pdf",width=14,height=10)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()
resSig = res[ res$padj < 0.1, ]
head( resSig[ order(resSig$pval), ] )
downreagulated
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )
up regulated 
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )
write.table(res, file="DE_seq_results.txt", sep="\t", quote=F, rownames=F)

## round values example
# seq1 <- c(1.256, 7.845, 7.735, 5.453, 1.945, 9.021, 4.265, 6.132, 7.549)
# mat1 <- matrix(seq1,3)
# mat1
#       [,1]  [,2]  [,3]
# [1,] 1.256 5.453 4.265
# [2,] 7.845 1.945 6.132
# [3,] 7.735 9.021 7.549
# round(mat1, 0)
#      [,1] [,2] [,3]
# [1,]    1    5    4
# [2,]    8    2    6
# [3,]    8    9    8
 
