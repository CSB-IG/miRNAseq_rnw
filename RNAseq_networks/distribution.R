# read squared matrix
w <- read.table(file = "x.adj.txt", header = T, sep = '\t')
colnames(w) <- rownames(w)
w <- data.matrix(frame = w, rownames.force = T)

# verify that is symmetric
isSymmetric(w)
q[upper.tri(q, diag=T)] <- NA

#graph quantiles
pdf("qq_plot.pdf",width=14,height=10)
qqnorm(w)
dev.off()

library(ggplot2)
library(reshape2)
x <- w
x[x == 0] <- NA
M1 <- log10(x)
ID <- as.data.frame(rownames(x))
M2 <- cbind(ID, M1)
l2 = melt(M2)
pdf("qq_plot.pdf",width=5,height=5)
ggplot(na.omit(l2), aes(value)) + geom_density(aes(group = variable), alpha=0.25)
dev.off()
