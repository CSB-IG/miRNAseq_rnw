##### FINAL SCRIPT ########### 

#filenames <- list.files("prueba", pattern="*.adj")
#set working dir = the dir where the .adj are

library(plyr)
library(minet)
library(igraph)

adjs <- list.files(pattern="*.adj")

#length must be equal to genes
length(adjs)

#read just the last line from file take values as characters,
#convert to matrix and assign first col as rownames
adjtosqmtx <- function(g){
  a <- scan(file = g, skip = 17, what = "")
  m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
  rownames(m) <- c(m[,1])
  m <- t(m[,2])
  m
}

#join all adjs in square matrix and fill empty spaces with NA
#replace all NA for 0
#order rownames and change character matrix to numeric
#select rownames as colnames as both are ordered alphabetically
M <- lapply(adjs, adjtosqmtx)
x <- t(rbind.fill.matrix(M))
x[is.na(x)] <- 0
x <- x[order(rownames(x)), ]
class(x) <- "numeric"
colnames(x) <- rownames(x)
#control: class must be 'matrix'
class(x)

#save squared matrix
write.table(x, file = "x.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

#convert to edgelist node-node-interaction (weighted=T) and save
xadj <- graph.adjacency(x,weighted=TRUE)
xsif <- get.data.frame(xadj)
write.table(xsif, file = "x.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

#apply data processing inequality (dpi)
xdpi <- aracne(mim =  x, eps = 0.2)

#save dpi squared matrix
write.table(xdpi, file = "xdpi.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

#convert dpi to edgelist node-node-interaction (weighted=T) and save
xdpiadj <- graph.adjacency(xdpi,weighted=TRUE)
xdpisif <- get.data.frame(xdpiadj)
write.table(xdpisif, file = "xdpi.sif.txt", sep = "\t", col.names= T, row.names= T, quote = F )

#load matrices as
z <- read.table(file = "matrixsqr_RNAseq752.txt", header = T, sep = '\t')


########
#repeated error aborted  
#obtain edgelist node-node-interaction
#library(PCIT)
#sifx <- getEdgeList(x, rm.zero=TRUE)
#sifxdpi <- getEdgeList(xdpi, rm.zero=TRUE)

