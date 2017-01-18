library(plyr)
library(minet)
library(igraph)

adjs <- sort(list.files(pattern="*_1.adj"))

length(adjs)

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

# sets colnames as the gene name from adjs name
names <- gsub('.{6}$','', adjs)
colnames(x) <- names

mir <- rownames(x)[grepl("hsa*", rownames(x))]
genes <- rownames(x)[grepl("^(?!hsa*)", rownames(x), perl=T)]

fullmtx <- function(x,mir,genes){
	r <- x[,mir]
	s <- x[mir,]
	if (sum(is.na(r)) < sum(is.na(s))){
		b <- rbind((x[genes,genes]), (t(r[genes,])))
		z <- cbind(r, (b[rownames(r),]))
		return(z)
		} else{
		b <- rbind((x[genes,genes]), (t(s[genes,])))
		z <- cbind(s, (b[rownames(s),]))
		return(z)}
}

x <- fullmtx(x,mir,genes)	

# replace all NA for 0
x[is.na(x)] <- 0

# turn character matrix to numeric
class(x) <- "numeric"

# select rownames as colnames and orders them alphabetically
x <- x[sort(rownames(x)), sort(colnames(x)) ]

isSymmetric(x)

# write.table(x, file = "basal_p1.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
# write.table(x, file = "her2_p1.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
# write.table(x, file = "luma_p1.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
# write.table(x, file = "lumb_p1.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
# write.table(x, file = "casos_p1.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
# write.table(x, file = "controles_p1.txt", sep = "\t", col.names= NA, row.names= T, quote = F )

# basal <- x en RData
# her2 <- x
# luma <- x
# lumb <- x
# control <- x
# caso <- x

basal <- data.matrix(read.table("basal_p1.txt"))
colnames(basal) <- rownames(basal)
her2 <- data.matrix(read.table("her2_p1.txt"))
colnames(her2) <- rownames(her2)
luma <- data.matrix(read.table("luma_p1.txt"))
colnames(luma) <- rownames(luma)
lumb <- data.matrix(read.table("lumb_p1.txt"))
colnames(lumb) <- rownames(lumb)
sanos <- data.matrix(read.table("controles_p1.txt"))
colnames(sanos) <- rownames(sanos)
casos <- data.matrix(read.table("casos_p1.txt"))
colnames(casos) <- rownames(casos)

Reduce(intersect, list(a,b,c))

nodos <- Reduce(intersect, list(colnames(luma), colnames(lumb), colnames(basal), colnames(her2), colnames(casos), colnames(sanos)))
# mir <- nodos[grepl("hsa*", nodos)]
# genes <- nodos[grepl("^(?!hsa*)", nodos, perl=T)]

basal <- basal[nodos, nodos]
her2 <- her2[nodos, nodos]
luma <- luma[nodos, nodos]
lumb <- lumb[nodos, nodos]
sanos <- sanos[nodos, nodos]
casos <- casos[nodos, nodos]

normalit<-function(m){
     (m - min(m[upper.tri(m)]))/(max(m)-min(m[upper.tri(m)]))
}

basal <- normalit(basal)
her2 <- normalit(her2)
luma <- normalit(luma)
lumb <- normalit(lumb)
casos <- normalit(casos)
sanos <- normalit(sanos)



write.table(basal, file = "basal_p1_filtro_norm.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
write.table(her2, file = "her2_p1_filtro_norm.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
write.table(luma, file = "luma_p1_filtro_norm.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
write.table(lumb, file = "lumb_p1_filtro_norm.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
write.table(casos, file = "casos_p1_filtro_norm.txt", sep = "\t", col.names= NA, row.names= T, quote = F )
write.table(sanos, file = "controles_p1_filtro_norm.txt", sep = "\t", col.names= NA, row.names= T, quote = F )


p2mi <- function(x, n){
  alfa <- 1.062
  beta <- -48.7
  gamma <- -0.634
  mi <- (alfa - log(x)) / ((-beta) + (-gamma)*n)
  return(mi)
  }












quantile(q, probs=seq(0.995, 1, 0.00001), na.rm=T)

