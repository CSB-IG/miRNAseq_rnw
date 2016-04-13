RNAseq <- read.table(file = "RNAseq_752_norm.rsem.aracne.txt", header=T, sep = '\t', row.names=1)
sanos <- read.table(file = "sanos.86.txt", header=T, sep = '\t', row.names=1)
length(intersect(rownames(RNAseq), rownames(sanos)))
RNAseq <- RNAseq[intersect(rownames(RNAseq), rownames(sanos)), ]
sanos <- sanos[intersect(rownames(RNAseq), rownames(sanos)), ]
identical(rownames(RNAseq), rownames(sanos))
identical(enfermos, sanos)
RNAseq <- cbind(RNAseq, sanos)
> library(plyr)
> library(NOISeq)
myfactors <- data.frame(Pacientes = c( rep('enfermos', 752), rep('sanos', 86)))
x <- readData(data= RNAseq, factors=myfactors)
QCreport(x, samples = NULL, factor = "Pacientes", norm = FALSE)


