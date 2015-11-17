#expression matrix
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

#turn NA to 0
#ejemplo: x[is.na(x)] <- 0
#M2[is.na(M2)] <- 0

#partir M2 en 10 pedazos
x1 <- M2[, 1:73]
x20 <- cbind(ID, x2)
x2 <- M2[, 74:147]
x20 <- cbind(ID, x2)
x3 <- M2[, 148:221]
x30 <- cbind(ID, x3)
x4 <- M2[, 222:295]
x40 <- cbind(ID, x4)
x5 <- M2[, 296:371]
x50 <- cbind(ID, x5)
x6 <- M2[, 372:445]
x60 <- cbind(ID, x6)
x7 <- M2[, 446:520]
x70 <- cbind(ID, x7)
x8 <- M2[, 521:593]
x80 <- cbind(ID, x8)
x9 <- M2[, 593:668]
x90 <- cbind(ID, x9)
x10 <- M2[, 669:753]
x100 <- cbind(ID, x10)

#input for density plots
lx1 = melt(x1)
lx2 = melt(x20)
lx3 = melt(x30)
lx4 = melt(x40)
lx5 = melt(x50)
lx6 = melt(x60)
lx7 = melt(x70)
lx8 = melt(x80)
lx9 = melt(x90)
lx10 = melt(x100)

#single plot
#ggplot(na.omit(l2), aes(x = value)) + geom_density(color = rainbow)

#overlayed plots
ggplot(na.omit(l2), aes(value)) + geom_density(aes(group = variable), alpha=0.25)

#different plots
ggplot(na.omit(lx10), aes (value)) + geom_density() + facet_wrap(~variable)

