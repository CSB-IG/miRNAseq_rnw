
library(plyr)

# scan tail of .adj to matrix  

a <- scan(file = "ZZZ3_1.adj", skip = 17, what = "")
m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
rownames(m) <- c(m[,1])
#dimnames(m) <- list(c(m[,1]), c("ID", a[1]))
m <- t(m[,2])

a1 <- scan(file = "DNAH2_1.adj", skip = 17, what = "")
m1 <- matrix(unlist(a1[-(1)]), ncol = 2, byrow = T)
dimnames(m1)<- list(c(m1[,1]), c("ID", a1[1]))
m1 <- t(m1[,2])

#bind matrix
x <- t(rbind.fill.matrix(m, m1))

# add colnames to matrix
colnames <- c(a[1], a1[1])

#sort matrix: mat <- mat[sort(rownames(mat)), sort(colnames(mat)) ]
x <- x[sort(rownames(x)), sort(colnames(x)) ]

# check if rownames in matrix and GENEID quadrate
# library(gplots)
# rowx <- rownames(x)
# listID <- scan("GENEID_filtro.txt", what = "", sep="\t")
# v <- venn(list(matrix_rownames = rowx,  ID_list = listID), show.plot = T)

# simplified example
# A <- matrix (1:4, 2)
# B <- matrix (6:11, 2)
# rbind.fill.matrix (A, B)

