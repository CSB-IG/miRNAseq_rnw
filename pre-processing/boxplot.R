
boxplot

bp <-cbind(T1, M[,3])
colnames(bp) <- c("ID", "variable", "TMM", "raw")

pdf("TMM-raw.miRNA.753.bp.pdf",width=10,height=10)
boxplot(bp[,3:4], data=bp)
dev.off()


setwd("/home/diana/")

enfermos <- read.table(file = "enfermos.86.txt", header = T, sep = '\t')
sanos <- read.table(file = "sanos.86.txt", header = T, sep = '\t')

boxplot(log((sanos)+1, 2))
boxplot(log((enfermos)+1, 2))
