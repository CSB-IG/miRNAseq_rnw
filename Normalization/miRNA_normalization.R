



library(limma)
library(edgeR)

# read miRNAs_752_maduros_filtro.txt
miRNA <- read.table("miRNAs_752_maduros_filtro.txt", sep="\t", header=T, row.names=1)

# x = expression matrix
# groups = factors(list)
# y = output

y = object: DGEList

object dispersion numeric vector of dispersion parameters.  By default, is extracted from
object or, if object contains no dispersion information, is set to 0.05
common.lib.size numeric scalar, the library size to normalize to; default is the geometric mean of
the original effective library sizes

> x <- read.delim("fileofcounts.txt",row.names="Symbol")
> group <- factor(c(1,1,2,2))
> y <- DGEList(counts=x,group=group)
> out <- equalizeLibSizes(object, dispersion=NULL, common.lib.size)
> y <- calcNormFactors(y, method="TMM, UQ, etc..")
> cpm.tmm <- cpm(y, normalized.lib.sizes=T)

"The effective library size is then the original library size multiplied by the scaling factor"

"common.lib.size numeric scalar, the library size to normalize to; default is the geometric mean of
the original effective library sizes"


# x = miRNA
group <- factor(c(rep('enfermos', 752)))
y <- DGEList(counts=miRNA,group=group)
# out <- equalizeLibSizes(y)
yTMM <- calcNormFactors(y, method="TMM")
yUQ <- calcNormFactors(y, method="upperquartile")

TMM <- cpm(yTMM, normalized.lib.sizes=T)
UQ <- cpm(yUQ, normalized.lib.sizes=T)

# TMM + UQ: x = TMM
TMM[is.na(TMM)] <- 0
y1TMM <- DGEList(counts=TMM,group=group)
y1TMMUQ <- calcNormFactors(y1TMM, method="upperquartile")  
TMMUQ <- cpm(y1TMMUQ, normalized.lib.sizes=T)

# UQ + TMM: x = UQ 
UQ[is.na(UQ)] <- 0
y1UQ <- DGEList(counts=UQ,group=group)
y1UQTMM <- calcNormFactors(y1UQ, method="TMM")  
UQTMM <- cpm(y1UQTMM, normalized.lib.sizes=T)

# density plot
library(ggplot2)
library(reshape2)

TMM <- as.data.frame(TMM)
UQ <- as.data.frame(UQ)
TMMUQ <- as.data.frame(TMMUQ)
UQTMM <- as.data.frame(UQTMM)

densityplot <- function(g){
  M1 <- log10(g+1)
  M2 <- cbind(rownames(M1), M1)
  rownames(M2) <- NULL
  colnames(M2)[1] <- "ID"
  l2 = melt(M2)
  l2
}


T1 <- densityplot(TMM)
U1 <- densityplot(UQ)
TU1 <- densityplot(TMMUQ)
UT1 <- densityplot(UQTMM)
M <- densityplot(miRNA)

pdf("TMM1.miRNA.753.pdf",width=10,height=10)
ggplot(T1, aes(value)) + geom_density(aes(group = variable))
dev.off()

pdf("UQ0.miRNA.753.pdf",width=10,height=10)
ggplot(U1, aes(value)) + geom_density(aes(group = variable))
dev.off()

pdf("TMMUQ0.miRNA.753.pdf",width=10,height=10)
ggplot(TU1, aes(value)) + geom_density(aes(group = variable))
dev.off()

pdf("UQTMM0.miRNA.753.pdf",width=10,height=10)
ggplot(UT1, aes(value)) + geom_density(aes(group = variable))
dev.off()
