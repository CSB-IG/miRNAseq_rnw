miRNAs_752_maduros_filtro.txt
library(ggplot2)
library(reshape2)
x <- read.table(file = "miRNAs_752_maduros_filtro.txt", header = T, sep = '\t', row.names=1)
> densityplot2 <- function(g){
+   M1 <- log2(g+1)
+   M2 <- cbind(rownames(M1), M1)
+   rownames(M2) <- NULL
+   colnames(M2)[1] <- "ID"
+   l2 = melt(M2)
+   l2
+ }

MIR752 <- densityplot2(x)

pdf("raw.miRNA.log2.752.pdf",width=10,height=10)
ggplot(MIR752, aes(value)) + geom_density(aes(group = variable))
dev.off()

xn <- read.table(file = "miRNA_752_normTMM.txt", header = T, sep = '\t', row.names=1)
MIR752N <- densityplot2(xn)
pdf("TMM.miRNA.log2.752.pdf",width=10,height=10)
ggplot(MIR752N, aes(value)) + geom_density(aes(group = variable))
dev.off()


me <- read.table(file = "miRNAs_enfermos_maduros_filtro.txt", header = T, sep = '\t', row.names=1)
ms <- read.table(file = "miRNAs_sanos_maduros_filtro.txt", header = T, sep = '\t', row.names=1)
identical(rownames(me), rownames(ms))
crudos <- cbind(me,ms)
enfermos <- densityplot2(me)
sanos <- densityplot2(ms)
mir <- densityplot2(crudos)
pdf("raw.miRNA.log2.enfermos.pdf",width=10,height=10)
ggplot(enfermos, aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("raw.miRNA.log2.sanos.pdf",width=10,height=10)
ggplot(sanos, aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("raw.miRNA.log2.86.pdf",width=10,height=10)
ggplot(mir, aes(value)) + geom_density(aes(group = variable))
dev.off()




mne <- read.table(file = "miRNAs_86_maduros_norm_enfermos.txt", header = T, sep = '\t', row.names=1)
mns <- read.table(file = "miRNAs_86_maduros_norm_sanos.txt", header = T, sep = '\t', row.names=1)
crudosn <- read.table(file = "miRNAs_86_maduros_norm.txt", header = T, sep = '\t', row.names=1)
enfermosn <- densityplot2(mne)
sanosn <- densityplot2(mns)
mirn <- densityplot2(crudosn)
pdf("TMM.miRNA.log2.enfermos.pdf",width=10,height=10)
ggplot(enfermosn, aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("TMM.miRNA.log2.sanos.pdf",width=10,height=10)
ggplot(sanosn, aes(value)) + geom_density(aes(group = variable))
dev.off()
pdf("TMM.miRNA.log2.86.pdf",width=10,height=10)
ggplot(mirn, aes(value)) + geom_density(aes(group = variable))
dev.off()


