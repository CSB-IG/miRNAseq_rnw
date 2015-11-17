# read data.frame
x <- read.table("data-frame.txt", sep = "\t")

# data.frame as numerical matrix
x <- data.matrix(frame = x, rownames.force = T)
is.matrix(x)
isSymmetrical(x)

# find elements that diverge between colnames and rownames

library(gplots)

a <- colnames(m)
b <- rownames(m)
c <- venn(list(columnas = a, celdas = b))

d <- cbind(a, b)
repeats <- data.frame(table(d))
repeats[which(repeats$Freq < 2 ),]

# change matrix colnames: colnames(matrix) <- c(lists)
colnames(m) <- c(b)

# compare gene list from condor & minet
w <- list.files(pattern="*.adj")
f <- venn(list(genes_condor = w,  genes_minet = b), show.plot = T)

# find genes that were ignored by condor but included in minet
setdiff(b, w)

# find genes that were included in condor but they are not genes
setdiff(w, b)

library(graph)
library(grid)
library(Rgraphviz)
library(igraph)

# convert squared matrix to graphNEL
gx <- as(x,"graphNEL")
gxdpi <- as(xdpi,"graphNEL")

# convert graphNEL to igraph
igx <- igraph.from.graphNEL(gx)
igxdpi <- igraph.from.graphNEL(gxdpi)

# convert squared matrix to igraph
xg <- graph.adjacency(adjmatrix= x, mode='undirected', diag=F, weighted=T)
xgdpi <- graph.adjacency(adjmatrix= xdpi, mode='undirected', diag=F, weighted=T)

# write as graphml
write.graph(igx, "matrix.graphml", format= "graphml")
write.graph(igxdpi, "matrix.graphml", format= "graphml")

