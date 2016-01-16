# read expression matrix
x <- read.csv("RNAseq_norm_753.csv")

library(ggplot2)
library(reshape2)

# turn 0 to NA: df[df == 0] <- NA
# x[x == 0] <- NA

# log10 transfomation and ID col
# M1 <- log10(x)

# log10 transfomation + 1 and ID col
# M1 <- log10(x)

# bind ID list from file 
# ID <- read.table("GeneID", header = TRUE)
# M2 <- cbind(ID, M1)

# bind ID list from data.frame 
# M2 <- cbind(rownames(M1), M1)

# turn NA to 0: x[is.na(x)] <- 0
# M2[is.na(M2)] <- 0

# nullily outlier: Matrix$col <- NULL
# M2$A245 <- NULL

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

x1 <- densityplot(x)

#single plot
#ggplot(na.omit(x1), aes(x = value)) + geom_density(color = rainbow)

#overlayed plots
ggplot(na.omit(x1), aes(value)) + geom_density(aes(group = variable), alpha=0.25)

#different plots
ggplot(na.omit(x1), aes (value)) + geom_density() + facet_wrap(~variable)

