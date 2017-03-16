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

# Get the location of the *.adj files
# Write the full path of the directory where the *.adj files are located
adj_files <- sort(list.files(file.choose(), full.names=TRUE, pattern="*.adj"))

# length must be equal to genes
# length(adjs)

# reads just the last line from file take values as characters
# converts to matrix and assign first column as rownames 

adjtosqmtx <- function(g){
  a <- scan(file = g, skip = 17, what = "")
  m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
  rownames(m) <- c(m[,1])
  m <- t(m[,2])
  m
}

# scale MI values from 0 to 1
normalit<-function(m){
     (m - min(m[upper.tri(m)]))/(max(m)-min(m[upper.tri(m)]))
}


build_mtx <- function(adj_files){
	#create matrix from individual adjs p-val=1
	M <- lapply(adj_files, adjtosqmtx)
	x <- t(rbind.fill.matrix(M))
	names <- gsub(".*/|_1.adj", "", adj_files)
	colnames(x) <- names
	# get mir and gene names
	mir <- rownames(x)[grepl("hsa*", rownames(x))]
	genes <- rownames(x)[grepl("^(?!hsa*)", rownames(x), perl=T)]
	#turn into a squared symmetric matrix
	# because MI networks are symmestrical the miR adjs were the only
	# ones that had the miR-gene edges, saving computing time
	r <- x[,mir]
	s <- x[mir,]
	if (sum(is.na(r)) < sum(is.na(s))){
		b <- rbind((x[genes,genes]), (t(r[genes,])))
		z <- cbind(r, (b[rownames(r),]))
		} else{
		b <- rbind((x[genes,genes]), (t(s[genes,])))
		z <- cbind(s, (b[rownames(s),]))
		}
	z[is.na(z)] <- 0
	class(z) <- "numeric"
	z <- z[sort(rownames(z)), sort(colnames(z))]
	z <- normalit(z)
	return(z)
}

x <- build_mtx(adj_files)

# save squared matrix
write.table(data.frame("Names"=rownames(x),x), file = "matrix_MI.txt", sep = "\t", col.names= FALSE, quote = F )

