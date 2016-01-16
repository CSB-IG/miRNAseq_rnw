# x <- expression matrix raw mature miRNAseq counts
# xn <- expression matrix normalized mature miRNAseq counts

library(ggplot2)
library(reshape2)

# turn 0 to NA: df[df == 0] <- NA
# x[x == 0] <- NA

#log10 transfomation and ID col
# M1 <- log10(x)

# set colnames as first column
# M2 <- cbind(rownames(M1), M1)

# delete rownames from data.frame
# rownames(M2) <- NULL

# turn NA to 0: x[is.na(x)] <- 0
# M2[is.na(M2)] <- 0

# nullify outlier: Matrix$col <- NULL
# M2$A245.01 <- NULL

# set variable name
# colnames(M2)[1] <- "ID"

# input for density plots
# l2 = melt(M2)

# density plot function
densityplot <- function(g){
  M1 <- log10(g+1)
  M2 <- cbind(rownames(M1), M1)
  rownames(M2) <- NULL
  colnames(M2)[1] <- "ID"
  l2 = melt(M2)
  l2
}

# x1 <- densityplot(x)

# single plot
# ggplot(na.omit(x1), aes(x = value)) + geom_density(color = rainbow)

#overlayed plots
pdf("raw.miRNA.753.pdf",width=10,height=10)
ggplot(na.omit(x1), aes(value)) + geom_density(aes(group = variable))
dev.off()

#different plots
#ggplot(na.omit(x1), aes (value)) + geom_density() + facet_wrap(~variable)
