library(plyr)
library(minet)
library(igraph)

adjs <- list.files(pattern="*_1.adj")
length(adjs)

adjtosqmtx <- function(g){
  a <- scan(file = g, skip = 17, what = "")
  m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
  rownames(m) <- c(m[,1])
  m <- t(m[,2])
  m
}

# 633 miRNAs
# enfermos



M <- lapply(adjs, adjtosqmtx)
x <- t(rbind.fill.matrix(M))
x[is.na(x)] <- 0
class(x) <- "numeric"
names <- gsub('.{6}$','', adjs)
colnames(x) <- names
x <- x[sort(rownames(x)), sort(colnames(x)) ]
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
x <- x[sort(rownames(x)), sort(colnames(x)) ]
class(x)
isSymmetric(x)
sum(x != 0)/2

# sum(x != 0)/2
# [1] 9734784

# enfermos
# 99.980%
# 0.14851600

p <- x
p[p < 0.14851600] <- 0
sum(p != 0)/2
# [1] 24867

pdpi <- aracne(mim =  p, eps = 0.1)
sum(p != 0)/2
dim(pdpi)
# sum(pdpi != 0)/2
# [1] 23582
# sum(diag(pdpi))
# [1] 0
# dim(pdpi)
[1] 15769 15769


write.table(pdpi, file = "enfermos.99.980.dpi.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

##### test if previous file had any errors
prueba <- read.table(file = "enfermos.mir.99.980.dpi.sif.txt", header = T, sep = '\t', stringsAsFactors=F)
xgdpi <- graph.adjacency(adjmatrix= pdpi, mode='undirected', diag=F, weighted=T)
test <- get.data.frame(xgdpi)
identical(prueba, test)
# TRUE

# xs es la matrix pequeña pero con el filtro de dpi y de p-value
xs <- pdpi[, names]
dim(xs)
# dim(xs)
# [1] 15769   633


enfermos <- read.table(file = "enfermos.99.87.corregido.dpi.txt", header = T, sep = '\t')

xs <- xs[sort(rownames(xs)), sort(colnames(xs)) ]
length(setdiff(rownames(xs), colnames(xs)))
genes <- setdiff(rownames(xs), colnames(xs))
y <- t(xs[genes,])
y <- y[ ,sort(colnames(y))]
enfermos <- enfermos[, sort(colnames(enfermos))]
colnames(enfermos) <- colnames(y)
identical(colnames(enfermos), colnames(y))
y <- rbind(y, enfermos)
y <- y[sort(rownames(y)),]
identical(rownames(xs), rownames(y))
x <-cbind(xs,y)
x <- x[sort(rownames(x)), sort(colnames(x)) ]
x <- data.matrix(x)
isSymmetric(x)
sum(pdpi != 0)/2
dim(pdpi)
sum(x != 0)/2
dim(x)

write.table(x, file = "enfermos.adj.completa.txt", sep = "\t", col.names= T, row.names= T, quote = F )


# sanos

matriz original###

adjs <- list.files(pattern="*_1.adj")
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
x <- x[sort(rownames(x)), sort(colnames(x)) ]
class(x)
x <- data.matrix(x)
isSymmetric(x)
sum(x != 0)/2

# sum(x != 0)/2
# [1] 9773668

# sanos
# 99.980%
# 0.2857994

p <- x
p[p < 0.2857994] <- 0

sum(p != 0)/2
# [1] 24865

pdpi <- aracne(mim =  p, eps = 0.1)
# sum(pdpi != 0)/2
# [1] 24401
# sum(diag(pdpi))
# [1] 0
# dim(pdpi)
# [1] 15769 15769
write.table(pdpi, file = "sanos.99.980.dpi.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

##### test if previous file had any errors
prueba <- read.table(file = "sanos.mir.99.980.dpi.sif.txt", header = T, sep = '\t', stringsAsFactors=F)
xgdpi <- graph.adjacency(adjmatrix= pdpi, mode='undirected', diag=F, weighted=T)
test <- get.data.frame(xgdpi)
identical(prueba, test)
# TRUE

# xs es la matrix pequeña pero con el filtro de dpi y de p-value
xs <- pdpi[, names]
dim(xs)
# dim(xs)
# [1] 15769   633

####

xs <- xs[sort(rownames(xs)), sort(colnames(xs)) ]


###### falta RNAseq sanos con dpi con nombre correcto
sanos <- read.table(file = "sanos.99.87.corregido.adj.txt", header = T, sep = '\t')
sanosdpi <- aracne(mim =  sanos, eps = 0.1)

dim(sanos)
dim(sanosdpi)
sum(sanos != 0)/2
sum(sanosdpi != 0)/2
write.table(sanosdpi, file = "sanos.99.987.dpi.rnaseq.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )
#####

length(colnames(pdpi))- length(colnames(sanosdpi))
# [1] 633

length(setdiff(rownames(xs), colnames(xs)))
genes <- setdiff(rownames(xs), colnames(xs))
y <- t(xs[genes,])
y <- y[ ,sort(colnames(y))]
sanosdpi <- sanosdpi[, sort(colnames(sanosdpi))]
colnames(sanosdpi) <- colnames(y)
identical(colnames(sanosdpi), colnames(y))
y <- rbind(y, sanosdpi)
y <- y[sort(rownames(y)),]
identical(rownames(xs), rownames(y))
x <-cbind(xs,y)
x <- x[sort(rownames(x)), sort(colnames(x)) ]
x <- data.matrix(x)
isSymmetric(x)
sum(pdpi != 0)/2
dim(pdpi)
sum(x != 0)/2
dim(x)

write.table(x, file = "sanos.adj.completa.txt", sep = "\t", col.names= T, row.names= T, quote = F )


























