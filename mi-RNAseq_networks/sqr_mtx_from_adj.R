# list of .adj to square matrix

library(plyr)
library(minet)
library(igraph)


# set working dir as the dir where the .adj are
# if other dir is used: filenames <- list.files("other_directory", pattern="*.adj")
adjs <- sort(list.files(pattern="*.adj"))

# length must be equal to genes
length(adjs)

# reads just the last line from file take values as characters
# converts to matrix and assign first col as rownames 

adjtosqmtx <- function(g){
  a <- scan(file = g, skip = 17, what = "")
  m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
  rownames(m) <- c(m[,1])
  m <- t(m[,2])
  m
}

# joins all adjs in square matrix and fill empty spaces with NA
M <- lapply(adjs, adjtosqmtx)
x <- t(rbind.fill.matrix(M))

# replace all NA for 0
x[is.na(x)] <- 0

# turn character matrix to numeric
class(x) <- "numeric"

# sets colnames as the gene name from adjs name
names <- gsub('.{6}$','', adjs)
colnames(x) <- names

# select rownames as colnames and orders them alphabetically
x <- x[sort(rownames(x)), sort(colnames(x)) ]

#### only for miRNA matrix or non squared matrix:
# must be equal to number of genes

length(setdiff(rownames(x), colnames(x)))

genes <- setdiff(rownames(x), colnames(x))
y <- t(x[genes,])
y2 <- mat.or.vec(length(colnames(y)), length(colnames(y)))
rownames(y2) <- colnames(y)
y <- rbind(y, y2)
x <- x[sort(rownames(x)),]
y <- y[sort(rownames(y)),]
identical(rownames(x), rownames(y))
x <-cbind(x,y)
####

# select rownames as colnames and orders them alphabetically
x <- x[sort(rownames(x)), sort(colnames(x)) ]

# control: class must be 'matrix'
# matrix must be symmetric
class(x)
isSymmetric(x)

# save squared matrix
write.table(x, file = "x.adj.txt", sep = "\t", col.names= NA, row.names= T, quote = F )

# convert to edgelist: node-node-interaction (weighted=T) and save
xg <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)
xsif <- get.data.frame(xg)
write.table(xsif, file = "x.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

# apply data processing inequality (dpi)
xdpi <- aracne(mim =  x, eps = 0.1)

#save dpi squared matrix
write.table(xdpi, file = "xdpi.adj.txt", sep = "\t", col.names= NA, row.names= T, quote = F )

#convert dpi to: edgelist node-node-interaction (weighted=T) and save
xgdpi <- graph.adjacency(adjmatrix= xdpi, mode='undirected', diag=F, weighted=T)
xdpisif <- get.data.frame(xgdpi)
write.table(xdpisif, file = "xdpi.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

#read matrices as:
x <- read.table(file = "x.adj.txt", header = T, sep = '\t')
x <- data.matrix(frame = x, rownames.force = T)

# whole script for source 

library(plyr)
library(minet)
library(igraph)

adjs <- sort(list.files(pattern="*.adj"))
length(adjs)

adjtosqmtx <- function(g){
  a <- scan(file = g, skip = 17, what = "")
  m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
  rownames(m) <- c(m[,1])
  m <- t(m[,2])
  m
}

M <- lapply(adjs, adjtosqmtx)
x <- t(rbind.fill.matrix(M))
x[is.na(x)] <- 0
class(x) <- "numeric"

names <- gsub('.{6}$','', adjs)
colnames(x) <- names
x <- x[sort(rownames(x)), sort(colnames(x)) ]

#### only for miRNA matrix or non squared matrix:
# must be equal to number of genes

length(setdiff(rownames(x), colnames(x)))

genes <- setdiff(rownames(x), colnames(x))
y <- t(x[genes,])
y2 <- mat.or.vec(length(colnames(y)), length(colnames(y)))
rownames(y2) <- colnames(y)
y <- rbind(y, y2)
x <- x[sort(rownames(x)),]
y <- y[sort(rownames(y)),]
identical(rownames(x), rownames(y))
x <-cbind(x,y)
####

x <- x[sort(rownames(x)), sort(colnames(x)) ]
class(x)
isSymmetric(x)
write.table(x, file = "x.adj.txt", sep = "\t", col.names= NA, row.names= T, quote = F )

xg <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)
xsif <- get.data.frame(xg)
write.table(xsif, file = "x.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

xdpi <- aracne(mim =  x, eps = 0.1)
write.table(xdpi, file = "xdpi.adj.txt", sep = "\t", col.names= NA, row.names= T, quote = F )

xgdpi <- graph.adjacency(adjmatrix= xdpi, mode='undirected', diag=F, weighted=T)
xdpisif <- get.data.frame(xgdpi)
write.table(xdpisif, file = "xdpi.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )
