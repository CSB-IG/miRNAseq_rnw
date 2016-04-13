> dimnames(x) <-list(c("a", "b"), c("c"))                                                                       
> x
  c
a 1
b 2

#matriz jueguete simetrica
> seq1 <- c(0, 1, 2, 1, 0, 1, 2, 1, 0)
> mat1 <- matrix(seq1,3)
> mat1
     [,1] [,2] [,3]
[1,]    0    1    2
[2,]    1    0    1
[3,]    2    1    0
> isSymmetric(mat1)
[1] TRUE


s<-matrix(1:25,5)
s[lower.tri(s)] = t(s)[lower.tri(s)]

 

seq1 <- seq(1:9)
mat1 <- matrix(seq1,3)
mat1
     [,1] [,2] [,3]
[1,]    1    4    7
[2,]    2    5    8
[3,]    3    6    9
> rownames(mat1) <- c('america', 'zeta', 'casa')
> colnames(mat1) <- rownames(mat1)
> mat1
      hola hola1 hola2
hola     1     4     7
hola1    2     5     8
hola2    3     6     9

#save matrix as
write.table(x, file = "matrixsqr_RNAseq752.txt", sep = "\t", col.names= T, row.names= T, quote = F )
write.table(y, file = "matrixsqr_RNAseq752.dpi.txt", sep = "\t", col.names= T, row.names= T, quote = F )

#load matrix as
z <- read.table(file = "matrixsqr_RNAseq752.dpi.txt", header = T, sep = '\t')
z <- as.matrix(z)

#filtrar por información mutua
mat1[mat1 < 0.1] <- 0

#####FINAL

library(plyr)
library(minet)
library(igraph)

adjs <- list.files(pattern="*.adj")
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
class(x)
isSymmetric(x)
write.table(x, file = "x.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xg <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)
xsif <- get.data.frame(xg)
write.table(xsif, file = "x.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

xdpi <- aracne(mim =  x, eps = 0.1)
write.table(xdpi, file = "xdpi.adj.txt", sep = "\t", col.names= T, row.names= F, quote = F )

xgdpi <- graph.adjacency(adjmatrix= xdpi, mode='undirected', diag=F, weighted=T)
xdpisif <- get.data.frame(xgdpi)
write.table(xdpisif, file = "xdpi.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

############

#### read matrix
x <- read.table(file = "x.adj.txt", header = T, sep = '\t')
x <- as.matrix(x)

tail -n +155437 x.sif.txt| awk '{print $2"\t"$3"\t"$4}' >> bueno.txt 


min val 0.00297737
max val 0.995974

awk '{print $2"\t"$3"\t"$4}' >> bueno.txt 

sifxdpi[which(sifxdpi$weight == 0.995974),]

01.nov.15

x <- read.table(file = "x.adj.txt", header = T, sep = '\t')
x <- as.matrix(x)
x[x < 0.1] <- 0
#x0.1 <- graph.adjacency(x,weighted=TRUE)
xg0.1 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)
xsif <- get.data.frame(xg)
write.table(xsif, file = "prueba1.txt", sep = "\t", col.names= T, row.names= F, quote = F)
x <- x[sort(rownames(x)), sort(colnames(x)) ]
a <- colnames(x)
rownames(x) <- a

isSymmetric(x)

03.nov.15
###########################
#hacer la que falta a p1
xgdpi <- graph.adjacency(adjmatrix= xdpi, mode='undirected', diag=F, weighted=T)
xdpisif <- get.data.frame(xgdpi)
write.table(xdpisif, file = "xdpi.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

#para terminar dpi a 0.2 de mi
write.table(xdpi02, file = "xdpi02.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )
xgdpi02 <- graph.adjacency(adjmatrix= xdpi02, mode='undirected', diag=F, weighted=T) 
xdpisif02 <- get.data.frame(xgdpi02)                                                                          
write.table(xdpisif02, file = "x.dpi.02.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )
#####################################3
# 0.3
x[x < 0.3] <- 0 
write.table(x, file = "x.03.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg03 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif03 <- get.data.frame(xg03)                                                                                                
write.table(xsif03, file = "x03.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi03 <- aracne(mim =  x, eps = 0.2)   
##### aqui me quede
write.table(xdpi03, file = "xdpi03.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi03 <- graph.adjacency(adjmatrix= xdpi03, mode='undirected', diag=F, weighted=T) 
xdpisif03 <- get.data.frame(xgdpi03)                                                                          
write.table(xdpisif03, file = "x.dpi.03.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

# 0.4
x[x < 0.4] <- 0 
write.table(x, file = "x.04.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg04 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif04 <- get.data.frame(xg04)                                                                                                
write.table(xsif04, file = "x04.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi04 <- aracne(mim =  x, eps = 0.2)   
write.table(xdpi04, file = "xdpi04.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi04 <- graph.adjacency(adjmatrix= xdpi04, mode='undirected', diag=F, weighted=T) 
xdpisif04 <- get.data.frame(xgdpi04)                                                                          
write.table(xdpisif04, file = "x.dpi.04.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

# 0.5
x[x < 0.5] <- 0 
write.table(x, file = "x.05.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg05 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif05 <- get.data.frame(xg05)                                                                                                
write.table(xsif05, file = "x05.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi05 <- aracne(mim =  x, eps = 0.2)   
write.table(xdpi05, file = "xdpi05.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi05 <- graph.adjacency(adjmatrix= xdpi05, mode='undirected', diag=F, weighted=T) 
xdpisif05 <- get.data.frame(xgdpi05)                                                                          
write.table(xdpisif05, file = "x.dpi.05.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )


# 0.6
x[x < 0.6] <- 0 
write.table(x, file = "x.06.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg06 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif06 <- get.data.frame(xg06)                                                                                                
write.table(xsif06, file = "x06.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi06 <- aracne(mim =  x, eps = 0.2)   
write.table(xdpi06, file = "xdpi06.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi06 <- graph.adjacency(adjmatrix= xdpi06, mode='undirected', diag=F, weighted=T) 
xdpisif06 <- get.data.frame(xgdpi06)                                                                          
write.table(xdpisif06, file = "x.dpi.06.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )


# 0.7
x[x < 0.7] <- 0 
write.table(x, file = "x.07.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg07 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif07 <- get.data.frame(xg07)                                                                                                
write.table(xsif07, file = "x07.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi07 <- aracne(mim =  x, eps = 0.2)   
write.table(xdpi07, file = "xdpi07.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi07 <- graph.adjacency(adjmatrix= xdpi07, mode='undirected', diag=F, weighted=T) 
xdpisif07 <- get.data.frame(xgdpi07)                                                                          
write.table(xdpisif07, file = "x.dpi.07.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )


# 0.8
x[x < 0.8] <- 0 
write.table(x, file = "x.08.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg08 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif08 <- get.data.frame(xg08)                                                                                                
write.table(xsif08, file = "x08.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi08 <- aracne(mim =  x, eps = 0.2)   
write.table(xdpi08, file = "xdpi08.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi08 <- graph.adjacency(adjmatrix= xdpi08, mode='undirected', diag=F, weighted=T) 
xdpisif08 <- get.data.frame(xgdpi08)                                                                          
write.table(xdpisif08, file = "x.dpi.08.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

# 0.9
x[x < 0.9] <- 0 
write.table(x, file = "x.09.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)
xg09 <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)                                  
xsif09 <- get.data.frame(xg09)                                                                                                
write.table(xsif09, file = "x09.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F ) 

xdpi09 <- aracne(mim =  x, eps = 0.2)   
write.table(xdpi09, file = "xdpi09.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F )

xgdpi09 <- graph.adjacency(adjmatrix= xdpi09, mode='undirected', diag=F, weighted=T) 
xdpisif09 <- get.data.frame(xgdpi09)                                                                          
write.table(xdpisif09, file = "x.dpi.09.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )

#repeated error aborted
#obtain edgelist node-node-interaction
#library(PCIT)
#sifx <- getEdgeList(x, rm.zero=TRUE)
#sifxdpi <- getEdgeList(xdpi, rm.zero=TRUE)

