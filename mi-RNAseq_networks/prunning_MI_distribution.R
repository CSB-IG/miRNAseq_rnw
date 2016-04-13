
# read p-value = 1 squared matrix
x <- read.table(file = "x.adj.txt", header = T, sep = '\t')
x <- data.matrix(frame = x, rownames.force = T)
colnames(x) <-rownames(x)

# x must be a numeric matrix and symmetric
class(x)
class(x[1,1])
isSymmetric(x)

# set a copy of x as q
q <- x

# set diagonal and upper triangle to NA
q[upper.tri(q, diag=T)] <- NA

# check that the upper triangle corresponds to NAs
q[1:10, 1:10]

# try different quantiles to find an appropiate number of interactions
quantile(q, na.rm=T)
quantile(q, probs=seq(0.9, 1, 0.001), na.rm=T)
quantile(q, probs=seq(0.995, 1, 0.00001), na.rm=T)

# change all the values smaller than the mi threshold to 0 in x
x[x < value] <- 0

# save filtered version of x
write.table(x, file = "x.quantile.adj.txt", sep = "\t", col.names= T, row.names= T, quote = F)

# save x as an edgelist: node node interaction

library(igraph)

xg <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)
xsif <- get.data.frame(xg)                                                                                    
write.table(xsif, file = "x.99.987.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )  
        
# apply different eps values to x at a determined mi cut off and save as an edgelist
                                              
library(minet)     
                                                                                           
xdpi <- aracne(mim =  x, eps = 0.15)                                                                          
xgdpi <- graph.adjacency(adjmatrix= xdpi, mode='undirected', diag=F, weighted=T)                              
xdpisif <- get.data.frame(xgdpi)                                                                              
write.table(xdpisif, file = "x.99.987.dpi.sif.txt", sep = "\t", col.names= T, row.names= F, quote = F )       



