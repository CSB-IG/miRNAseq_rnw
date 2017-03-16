################################################################################
# Normalize miRNAs by TMM using the edgeR package (Robinson et at. 2010)
# Diana Drago García's - Github: https://github.com/CSB-IG/miRNAseq_rnw
################################################################################

# Install/load package ##############

library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)

################################################
# Load expression data for miRNA normalization #
# The miRNA file should have the following     #
# format:                                      #
#                                              #
# NAME	Sample1	Sample2	Sample3	...            #
# miRNA1	value	value	value	...        #
# miRNA2	value	value	value	...        #
################################################

# miRNA_data = expression matrix
# groups = factors(list)
# y = output
# the output is an object of the DGEList type. y = object: DGEList

# read miRNA data file
miRNA_data <- read.table("miRNA_data.txt", sep="\t", header=T, row.names=1)

# The numbers should be replaced with the number of tumour and control samples
# according to the samples in the miRNA_data matrix
# eg. sample1 = tumour, sample2=tumour, sample3=control,...
group <- factor(c("tumour","tumour","control"...))


miRNA_DGEList <- DGEList(counts=miRNA_data,group=group)

# TMM normalization with default parameters
TMM_normFact <- calcNormFactors(miRNA_DGEList, method="TMM")

# CPM expression matrix
TMM <- cpm(TMM_normFact, normalized.lib.sizes=T)

# density plot

TMM <- as.data.frame(TMM)

densityplot2 <- function(g){
+   M1 <- log2(g+1)
+   M2 <- cbind(rownames(M1), M1)
+   rownames(M2) <- NULL
+   colnames(M2)[1] <- "ID"
+   l2 = melt(M2)
+   l2
+ }

TMM_density <- densityplot(TMM)

pdf("TMM_miRNA.pdf",width=10,height=10)
ggplot(TMM_density, aes(value)) + geom_density(aes(group = variable))
dev.off()
