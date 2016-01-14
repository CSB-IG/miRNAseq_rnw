library(gplots)
library(RColorBrewer)

genes <- scan("DE_genes_Deseq2.txt", what="")
dataset <- read.table("RNAseq_752_norm.rsem.aracne.txt", header=T, sep="\t", row.names=1)
subtypes <- read.table("Subtipos_752_RNAseq.txt", header=T, sep="\t", row.names=1)
dataset <- data.matrix(dataset)
dataset[dataset == 0] <- NA
dataset <- log10(dataset)
dataset[is.na(dataset)] <- 0

rows.to.keep<-which(rownames(dataset) %in% genes)
dataset <- dataset[rows.to.keep,]


pacientes <- as.vector(row.names(subtypes))
dataset <- dataset[,match(pacientes, colnames(dataset))]

orden <- subtypes[,1]
orden <- gsub("LumA", "darkblue", orden)
orden <- gsub("LumB", "firebrick2", orden)
orden <- gsub("Her2", "orange", orden)
orden <- gsub("Basal", "yellow", orden)
orden <- gsub("Normal", "steelblue1", orden)


my_palette <- colorRampPalette(c("royalblue", "linen", "red3"))(n = 299)

# for un even breaks
col_breaks = c(seq(-2.366532,1.0,length=100), seq(1.1,3.5,length=100), seq(3.6,6.264798,length=100))

# squared matrix
distance = dist(dataset, method = "euclidean")
cluster = hclust(distance, method = "complete")

# not squared
row_distance = dist(dataset, method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
col_distance = dist(t(dataset), method = "euclidean")
col_cluster = hclust(col_distance, method = "complete")


#pdf("heatmap1.pdf",width=3,height=3)
#png("heatmap1.png", width = 6*500, height = 6*400, res = 300, pointsize = 10)
pdf("heatmap1.pdf",width=5,height=5)     
heatmap.2(x = dataset,  
  main = "RNAseq",
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins = c(12,9),     
  col = my_palette, 
  breaks = col_breaks,
  ColSideColors = orden,          
  Rowv = as.dendrogram(row_cluster), 
  Colv = as.dendrogram(col_cluster))       
dev.off()


ColSideColors = c(   
     rep("gray", 3),   # categories, Measurement 1-3: green
     rep("blue", 3),    # Measurement 4-6: blue
     rep("black", 4)),    # Measurement 7-10: red


heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
  breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA")            # turn off column clustering



RowSideColors = c(    # grouping row-variables into different
     rep("gray", 3),   # categories, Measurement 1-3: green
     rep("blue", 3),    # Measurement 4-6: blue
     rep("black", 4)),    # Measurement 7-10: red

######
dataset <- read.table("DE_seq2_results.txt", header=T, sep="\t", row.names=1)
res <- dataset
sum(res$padj < 0.05, na.rm=TRUE)
arriba <- genes1[genes1$log2FoldChange > 1.5,]
abajo <- genes1[genes1$log2FoldChange < -1.5,]
lista <- rbind(arriba, abajo)
write.table(lista, file="DE_seq2_results_dif_exp_pvalue.txt", sep="\t", quote=F)

