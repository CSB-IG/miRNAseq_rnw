# read squared matrix
w <- read.table(file = "x.adj.txt", header = T, sep = '\t')
colnames(w) <- rownames(w)
w <- data.matrix(frame = w, rownames.force = T)

# verify that is symmetric
isSymmetric(w)
q[upper.tri(q, diag=T)] <- NA

#graph quantiles
qqnorm(w)
