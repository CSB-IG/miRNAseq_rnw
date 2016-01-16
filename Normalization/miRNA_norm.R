



library(limma)
library(edgeR)
read miRNAs_752_maduros_filtro.txt
miRNA <- read.table("miRNAs_752_maduros_filtro.txt", sep="\t", header=T, row.names=1)

# x = expression matrix
# groups = factors(list)
# y = output

equal library sizes
object: DGEList
object dispersion numeric vector of dispersion parameters.  By default, is extracted from
object or, if object contains no dispersion information, is set to 0.05
common.lib.size numeric scalar, the library size to normalize to; default is the geometric mean of
the original effective library sizes

> x <- read.delim("fileofcounts.txt",row.names="Symbol")
> group <- factor(c(1,1,2,2))
> y <- DGEList(counts=x,group=group)
> out <- equalizeLibSizes(object, dispersion=NULL, common.lib.size)
> y <- calcNormFactors(y)

"The effective library size is then the original library size multiplied by the scaling factor"

"common.lib.size numeric scalar, the library size to normalize to; default is the geometric mean of
the original effective library sizes"

normcountsTMM <- 1e06*t(t(d$counts) / (d$samples$lib.size*d$samples$norm.factors) )

normcountsTMM <- 1e06*((miRNA / (yTMM$samples$lib.size*yTMM$samples$norm.factors) )


# x = miRNA
group <- factor(c(rep('enfermos', 752)))
y <- DGEList(counts=miRNA,group=group)
out <- equalizeLibSizes(y)
yTMM <- calcNormFactors(y, method="TMM")
yUQ <- calcNormFactors(y, method="upperquartile")

# normfactors*matrix: t(t(m) * v)
TMM <- t(t(miRNA) * yTMM$samples$norm.factors)
UQ <- t(t(miRNA) * yUQ$samples$norm.factors)

# TMM + UQ: x = TMM
TMM[is.na(TMM)] <- 0
y1TMM <- DGEList(counts=TMM,group=group)
y1TMMUQ <- calcNormFactors(y1TMM, method="upperquartile")  
TMMUQ <- t(t(TMM) * y1TMMUQ$samples$norm.factors)

# UQ + TMM: x = UQ 
UQ[is.na(UQ)] <- 0
y1UQ <- DGEList(counts=UQ,group=group)
y1UQTMM <- calcNormFactors(y1UQ, method="upperquartile")  
UQTMM <- t(t(UQ) * y1UQTMM$samples$norm.factors)

library(ggplot2)
library(reshape2)

TMM <- as.data.frame(TMM)
UQ <- as.data.frame(UQ)
TMMUQ <- as.data.frame(TMMUQ)
UQTMM <- as.data.frame(UQTMM)

TMM[TMM == 0] <- NA
UQ[UQ == 0] <- NA
TMMUQ[TMMUQ == 0] <- NA
UQTMM[UQTMM == 0] <- NA
TM1 <- log10(TMM)
UM1 <- log10(UQ)
TUM1 <- log10(TMMUQ)
UTM1 <- log10(UQTMM)
TM2 <- cbind(rownames(TM1), TM1)
UM2 <- cbind(rownames(UM1), UM1)
TUM2 <- cbind(rownames(TUM1), TUM1)
UTM2 <- cbind(rownames(UTM1), UTM1)
rownames(TM2) <- NULL
rownames(UM2) <- NULL
rownames(TUM2) <- NULL
rownames(UTM2) <- NULL
colnames(TM2)[1] <- "ID"
colnames(UM2)[1] <- "ID"
colnames(TUM2)[1] <- "ID"
colnames(UTM2)[1] <- "ID"
Tl2 = melt(TM2)
Ul2 = melt(UM2)
TUl2 = melt(TUM2)
UTl2 = melt(UTM2)
pdf("TMM.miRNA.753.pdf",width=10,height=10)
ggplot(na.omit(Tl2), aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("UQ.miRNA.753.pdf",width=10,height=10)
ggplot(na.omit(Ul2), aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("TMMUQ.miRNA.753.pdf",width=10,height=10)
ggplot(na.omit(TUl2), aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("UQTMM.miRNA.753.pdf",width=10,height=10)
ggplot(na.omit(UTl2), aes(value)) + geom_density(aes(group = variable))
dev.off()
final regresar NA a ceros!!!!!