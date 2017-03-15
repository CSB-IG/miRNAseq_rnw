# from list of .adj to a square matrix, graph and edgelist 

################################################
# Directory with all files *_1.adj             #
# *_1.adj files have the format:               #
#                                              #
# >  Input file      /input_file.txt           #
# >  ADJ file                                  #
# >  Output file     /output_file_1.adj        #
# >  Algorithm       fixed_bandwidth           #
# >  Kernel width    0.180679                  #
# >  MI threshold    0                         #
# >  MI P-value      1                         #
# >  DPI tolerance   1                         #
# >  Correction      0                         #
# >  Subnetwork file                           #
# >  Hub probe       _Gene1                    #
# >  Control probe                             #
# >  Condition                                 #
# >  Percentage      0.35                      #
# >  TF annotation                             #
# >  Filter mean     0                         #
# >  Filter CV       0                         #
# Gene1	Gene2	value	Gene3	value	...    #
################################################

# Install/load packages ##############

library(plyr)
library(minet)
library(igraph)

# set the directory where the *.adj files are located as the working directory
# if other dir is used: filenames <- list.files("other_directory", pattern="*.adj")
adjs <- sort(list.files(pattern="*_1.adj"))

# length must be equal to genes
length(adjs)

# reads just the last line from file take values as characters
# converts to matrix and assign first column as rownames 

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

# sets colnames as the gene/miR name from adjs file name in the directory
names <- gsub('.{6}$','', adjs)
colnames(x) <- names

# get miR and gene name lists
mir <- rownames(x)[grepl("hsa*", rownames(x))]
genes <- rownames(x)[grepl("^(?!hsa*)", rownames(x), perl=T)]

# makes matrices symmetric with the triangle that has more information
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

# get squared matrix
x <- fullmtx(x,mir,genes)	

# replace all NA for 0
x[is.na(x)] <- 0

# turn character matrix to numeric
class(x) <- "numeric"

# select rownames as colnames and orders them alphabetically
x <- x[sort(rownames(x)), sort(colnames(x)) ]

isSymmetric(x)

# scale MI values from 0 to 1
normalit<-function(m){
     (m - min(m[upper.tri(m)]))/(max(m)-min(m[upper.tri(m)]))
}

x <- normalit(x)

# save squared matrix
write.table(x("Names"=rownames(x),x), file = "matrix_MI.txt", sep = "\t", col.names= NA, row.names= T, quote = F )

