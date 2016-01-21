
boxplot

bp <-cbind(T1, M[,3])
colnames(bp) <- c("ID", "variable", "TMM", "raw")

pdf("TMM-raw.miRNA.753.bp.pdf",width=10,height=10)
boxplot(bp[,3:4], data=bp)
dev.off()