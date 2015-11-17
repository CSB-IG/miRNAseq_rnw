# read expression matrix
x <- read.csv("RNAseq_norm_753.csv")

library(ggplot2)
library(reshape2)

#turn 0 to NA: df[df == 0] <- NA
x[x == 0] <- NA

#log10 transfomation and ID col
M1 <- log10(x)

# bind ID list
ID <- read.table("GeneID", header = TRUE)
M2 <- cbind(ID, M1)

#turn NA to 0: x[is.na(x)] <- 0
#M2[is.na(M2)] <- 0

#nullily outlier: Matrix$col <- NULL
M2$A245 <- NULL

#input for density plots
l2 = melt(M2)

#single plot
#ggplot(na.omit(l2), aes(x = value)) + geom_density(color = rainbow)

#overlayed plots
ggplot(na.omit(l2), aes(value)) + geom_density(aes(group = variable), alpha=0.25)

#different plots
ggplot(na.omit(l2), aes (value)) + geom_density() + facet_wrap(~variable)

