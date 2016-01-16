library(minet)
# applies the data processing inequality to triplets of nodes

x <- read.table("squared_matrix.txt", header=T, sep="\t")
x <- as.matrix(x)
isSymmetric(x)

ptm <- proc.time()
xdpi <- aracne(x, eps=0.1)
proc.time() - ptm

# returns weighted adjacency matrix of the network

dataset <- read.table("expression_matrix.txt", header=T, sep="\t")
dataset <- as.matrix(dataset)
dataset <- t(dataset)
dataset <- as.data.frame(dataset)

ptm <- proc.time()
RNAseq_aracne <- minet(dataset, method="aracne", estimator="mi.empirical", disc="globalequalwidth")
proc.time() - ptm

# build mutual information matrix

ptm <- proc.time()
mim <- build.mim(dataset,estimator="mi.empirical", disc = "globalequalwidth")
dpi <- aracne(mim, eps = 0.15)
proc.time() - ptm
