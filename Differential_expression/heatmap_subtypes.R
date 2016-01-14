library(gplots)
library(RColorBrewer)

#genes DE
basal <- read.table("DE_basal_control.txt", header=T, sep="\t", row.names=1)
her2 <- read.table("DE_her2_control.txt", header=T, sep="\t", row.names=1)
luma <- read.table("DE_luma_control.txt", header=T, sep="\t", row.names=1)
lumb <- read.table("DE_lumb_control.txt", header=T, sep="\t", row.names=1)
normal <- read.table("DE_normal_control.txt", header=T, sep="\t", row.names=1)

basal_g <- rbind(basal[which(basal$log2FoldChange > 1.5 & basal$padj < 0.05),], 
		basal[which(basal$log2FoldChange < -1.5 & basal$padj < 0.05),])

her2_g <- rbind(her2[which(her2$log2FoldChange > 1.5 & her2$padj < 0.05),], 
		her2[which(her2$log2FoldChange < -1.5 & her2$padj < 0.05),])

luma_g <- rbind(luma[which(luma$log2FoldChange > 1.5 & luma$padj < 0.05),], 
		luma[which(luma$log2FoldChange < -1.5 & luma$padj < 0.05),])

lumb_g <- rbind(lumb[which(lumb$log2FoldChange > 1.5 & lumb$padj < 0.05),], 
		lumb[which(lumb$log2FoldChange < -1.5 & lumb$padj < 0.05),])

normal_g <- rbind(normal[which(normal$log2FoldChange > 1.5 & normal$padj < 0.05),], 
		normal[which(normal$log2FoldChange < -1.5 & normal$padj < 0.05),])

genes <-Reduce(union, list(rownames(basal_g), rownames(her2_g), rownames(luma_g), rownames(lumb_g), rownames(normal_g)))

dataset <- read.table("RNAseq_752_norm.rsem.aracne.txt", header=T, sep="\t", row.names=1)
subtypes <- read.table("Subtipos_pam50.scale_752_RNAseq.txt", sep="\t", row.names=1, stringsAsFactors=F)
dataset <- data.matrix(dataset)
dataset[dataset == 0] <- NA
dataset <- log10(dataset)
dataset[is.na(dataset)] <- 0

rows.to.keep <-which(rownames(dataset) %in% genes)
dataset <- dataset[rows.to.keep,]

subtypes <- subtypes[ order(subtypes[,1]), , drop=FALSE]

pacientes <- as.vector(row.names(subtypes))
dataset <- dataset[,match(pacientes, colnames(dataset))]

orden <- subtypes[,1]
orden <- gsub("LumA", "darkblue", orden)
orden <- gsub("LumB", "firebrick2", orden)
orden <- gsub("Her2", "orange", orden)
orden <- gsub("Basal", "yellow", orden)
orden <- gsub("Normal", "steelblue1", orden)


my_palette <- colorRampPalette(c("royalblue", "linen", "red3"))(n = 299)

min() max()
# for un even breaks
col_breaks = c(seq(-2.39794,1.0,length=100), seq(1.1,3.5,length=100), seq(3.6,6.315139,length=100))

# squared matrix
distance = dist(dataset, method = "euclidean")
cluster = hclust(distance, method = "complete")

# not squared
row_distance = dist(dataset, method = "manhattan")
row_cluster = hclust(row_distance, method = "median")
col_distance = dist(t(dataset), method = "manhattan")
col_cluster = hclust(col_distance, method = "median")


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


pdf("heatmap3.pdf",width=5,height=5)     
heatmap.2(x = dataset,  
  main = "RNAseq",
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins = c(12,9),     
  col = my_palette, 
  breaks = col_breaks,
  ColSideColors = orden,
  labRow=F,
  labCol=F)                  
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


#####################




basal_up <- rownames(basal[which(basal$log2FoldChange > 1.5 & basal$padj < 0.05),])
basal_down <- rownames(basal[which(basal$log2FoldChange < -1.5 & basal$padj < 0.05),])
her2_up <- rownames(her2[which(her2$log2FoldChange > 1.5 & her2$padj < 0.05),])		
her2_down <- rownames(her2[which(her2$log2FoldChange < -1.5 & her2$padj < 0.05),])
luma_up <- rownames(luma[which(luma$log2FoldChange > 1.5 & luma$padj < 0.05),]) 	
luma_down <- rownames(luma[which(luma$log2FoldChange < -1.5 & luma$padj < 0.05),])
lumb_up <- rownames(lumb[which(lumb$log2FoldChange > 1.5 & lumb$padj < 0.05),]) 		
lumb_down <- rownames(lumb[which(lumb$log2FoldChange < -1.5 & lumb$padj < 0.05),])
normal_up <- rownames(normal[which(normal$log2FoldChange > 1.5 & normal$padj < 0.05),]) 	
normal_down <- rownames(normal[which(normal$log2FoldChange < -1.5 & normal$padj < 0.05),])


basal_up
her2_up	
luma_up	
lumb_up		
normal_up


basal_down
her2_down
luma_down
lumb_down
normal_down

pdf("venn_basal_up.pdf",width=5,height=5)
venn(list(basal_up=basal_up, her2_down=her2_down, luma_down=luma_down, lumb_down=lumb_down, normal_down=normal_down))
dev.off()

b <- Reduce(intersect, list(basal_up, her2_down, luma_down, lumb_down, normal_down))

pdf("venn_her2_up.pdf",width=5,height=5)
venn(list(her2_up=her2_up, basal_down=basal_down, luma_down=luma_down, lumb_down=lumb_down, normal_down=normal_down))
dev.off()

pdf("venn_luma_up.pdf",width=5,height=5)
venn(list(luma_up=luma_up, basal_down=basal_down, her2_down=her2_down, lumb_down=lumb_down, normal_down=normal_down))
dev.off()

c <- Reduce(intersect, list(luma_up, her2_down, basal_down, lumb_down, normal_down))

pdf("venn_lumb_up.pdf",width=5,height=5)
venn(list(lumb_up=lumb_up, basal_down=basal_down, her2_down=her2_down, luma_down=luma_down, normal_down=normal_down))
dev.off()

pdf("venn_normal_up.pdf",width=5,height=5)
venn(list(normal_up=normal_up, basal_down=basal_down, her2_down=her2_down, luma_down=luma_down, lumb_down=lumb_down))
dev.off()




pdf("venn_basal_down.pdf",width=5,height=5)
venn(list(basal_down=basal_down, her2_up=her2_up, luma_up=luma_up, lumb_up=lumb_up, normal_up=normal_up))
dev.off()

a <- Reduce(intersect, list(basal_down, her2_up, luma_up, lumb_up, normal_up))

pdf("venn_her2_down.pdf",width=5,height=5)
venn(list(her2_down=her2_down, basal_up=basal_up, luma_up=luma_up, lumb_up=lumb_up, normal_up=normal_up))
dev.off()

pdf("venn_luma_down.pdf",width=5,height=5)
venn(list(luma_down=luma_down, basal_up=basal_up, her2_up=her2_up, lumb_up=lumb_up, normal_up=normal_up))
dev.off()

pdf("venn_lumb_down.pdf",width=5,height=5)
venn(list(lumb_down=lumb_down, basal_up=basal_up, her2_up=her2_up, luma_up=luma_up, normal_up=normal_up))
dev.off()

pdf("venn_normal_down.pdf",width=5,height=5)
venn(list(normal_down=normal_down, basal_up=basal_up, her2_up=her2_up, luma_up=luma_up, lumb_up=lumb_up))
dev.off()


d <- c(a, b, c)

rows.to.keep2 <-which(rownames(dataset) %in% d)
dataset2 <- dataset[rows.to.keep2,]


subtypes <- subtypes[ order(subtypes[,1]), , drop=FALSE]

pacientes <- as.vector(row.names(subtypes))
dataset2 <- dataset2[,match(pacientes, colnames(dataset2))]

orden <- subtypes[,1]
orden <- gsub("LumA", "darkblue", orden)
orden <- gsub("LumB", "firebrick2", orden)
orden <- gsub("Her2", "orange", orden)
orden <- gsub("Basal", "yellow", orden)
orden <- gsub("Normal", "steelblue1", orden)


my_palette <- colorRampPalette(c("royalblue", "linen", "red3"))(n = 299)

min() max()
# for un even breaks
#col_breaks2 = c(seq(-0.7311881,0.0,length=100), seq(0.001,2.4,length=100), seq(2.41, 5.66473,length=100))
#col_breaks2 = c(seq(-0.7311881,0.0,length=100), seq(0.01,0.56,length=100), seq(0.57, 5.66473,length=100)) 
# squared matrix
distance = dist(dataset, method = "euclidean")
cluster = hclust(distance, method = "complete")

# not squared
row_distance2 = dist(dataset2, method = "manhattan")
row_cluster2 = hclust(row_distance2, method = "median")
col_distance2 = dist(t(dataset2), method = "manhattan")
col_cluster2 = hclust(col_distance2, method = "median")

col_breaks2 = c(seq(-0.7311881,-0.451,length=100), seq(-0.45,2.2,length=100), seq(2.21, 5.66473,length=100))
col_breaks2 = c(seq(-0.7311881,-0.491,length=100), seq(-0.49,1.3,length=100), seq(1.35, 5.66473,length=100))

#pdf("heatmap1.pdf",width=3,height=3)
#png("heatmap1.png", width = 6*500, height = 6*400, res = 300, pointsize = 10)
pdf("heatmap4.pdf",width=8,height=8)     
heatmap.2(x = dataset2,  
  main = "RNAseq",
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins = c(12,9),     
  col = my_palette, 
  breaks = col_breaks2,
  ColSideColors = orden, 
  Colv=FALSE,   
  labCol=FALSE,      
  Rowv = as.dendrogram(row_cluster2)) 
 
par(lend = 1)           
legend("topright",      
       legend = c("LumA", "LumB", "Her2", "Basal", "Normal"), 
       col = c("darkblue", "firebrick2", "orange", "yellow", "steelblue"), 
       lty= 1,         
       lwd = 5          
)

    
dev.off()


################# sin normal like


pdf("venn_basal_up1.pdf",width=10,height=7)
venn(list(basal_up=basal_up, her2_down=her2_down, luma_down=luma_down, lumb_down=lumb_down))
dev.off()

b1 <- Reduce(intersect, list(basal_up, her2_down, luma_down, lumb_down))

pdf("venn_her2_up1.pdf",width=10,height=7)
venn(list(her2_up=her2_up, basal_down=basal_down, luma_down=luma_down, lumb_down=lumb_down))
dev.off()

pdf("venn_luma_up1.pdf",width=10,height=7)
venn(list(luma_up=luma_up, basal_down=basal_down, her2_down=her2_down, lumb_down=lumb_down))
dev.off()

c1 <- Reduce(intersect, list(luma_up, her2_down, basal_down, lumb_down))

pdf("venn_lumb_up1.pdf",width=10,height=7)
venn(list(lumb_up=lumb_up, basal_down=basal_down, her2_down=her2_down, luma_down=luma_down))
dev.off()




pdf("venn_basal_down1.pdf",width=10,height=7)
venn(list(basal_down=basal_down, her2_up=her2_up, luma_up=luma_up, lumb_up=lumb_up))
dev.off()

a1 <- Reduce(intersect, list(basal_down, her2_up, luma_up, lumb_up))

pdf("venn_her2_down1.pdf",width=10,height=7)
venn(list(her2_down=her2_down, basal_up=basal_up, luma_up=luma_up, lumb_up=lumb_up))
dev.off()

pdf("venn_luma_down1.pdf",width=10,height=7)
venn(list(luma_down=luma_down, basal_up=basal_up, her2_up=her2_up, lumb_up=lumb_up))
dev.off()

pdf("venn_lumb_down1.pdf",width=10,height=7)
venn(list(lumb_down=lumb_down, basal_up=basal_up, her2_up=her2_up, luma_up=luma_up))
dev.off()


d1 <- c(a1, b1, c1)

rows.to.keep3 <-which(rownames(dataset) %in% d1)
dataset3 <- dataset[rows.to.keep3,]


subtypes <- subtypes[ order(subtypes[,1]), , drop=FALSE]

pacientes <- as.vector(row.names(subtypes))
pacientes2 <- pacientes[1:688]
dataset3 <- dataset3[,match(pacientes2, colnames(dataset3))]

orden <- subtypes[,1]
orden <- gsub("LumA", "darkblue", orden)
orden <- gsub("LumB", "firebrick2", orden)
orden <- gsub("Her2", "orange", orden)
orden <- gsub("Basal", "yellow", orden)
orden <- gsub("Normal", "steelblue1", orden)
orden2 <- orden[1:688]

my_palette <- colorRampPalette(c("royalblue", "linen", "red3"))(n = 299)

min() max()
# for un even breaks
col_breaks3 = c(seq(-0.7311881,0.0,length=100), seq(0.01,0.26,length=100), seq(0.261, 5.919181,length=100)) 

col_breaks3 = c(seq(-0.7311881,0.1,length=100), seq(0.11,0.36,length=100), seq(0.361, 5.919181,length=100)) 
col_breaks3 = c(seq(-0.7311881,0.15,length=100), seq(0.151,0.36,length=150), seq(0.361, 5.919181,length=50))

col_breaks3 = c(seq(-0.7311881,0.24,length=100), seq(0.25, 1.5,length=100), seq(1.51, 5.919181,length=100))
# not squared
row_distance3 = dist(dataset3, method = "manhattan")
row_cluster3 = hclust(row_distance3, method = "median")
col_distance3 = dist(t(dataset3), method = "manhattan")
col_cluster3 = hclust(col_distance3, method = "median")

#pdf("heatmap1.pdf",width=3,height=3)
#png("heatmap1.png", width = 6*500, height = 6*400, res = 300, pointsize = 10)
pdf("heatmap52.pdf",width=8,height=8)     
heatmap.2(x = dataset3,  
  main = "RNAseq",
  notecol="black",      
  density.info="none",  
  trace="none",         
  margins = c(12,9),     
  col = my_palette, 
  breaks = col_breaks3,
  ColSideColors = orden2, 
  Colv=FALSE,   
  labCol=FALSE,      
  Rowv = as.dendrogram(row_cluster3)) 
 
par(lend = 1)           
legend("topright",      
       legend = c("LumA", "LumB", "Her2", "Basal"), 
       col = c("darkblue", "firebrick2", "orange", "yellow"), 
       lty= 1,         
       lwd = 5          
)

    
dev.off()
 #### con log fold

basal <- read.table("DE_basal_control.txt", header=T, sep="\t", row.names=1)
her2 <- read.table("DE_her2_control.txt", header=T, sep="\t", row.names=1)
luma <- read.table("DE_luma_control.txt", header=T, sep="\t", row.names=1)
lumb <- read.table("DE_lumb_control.txt", header=T, sep="\t", row.names=1)
normal <- read.table("DE_normal_control.txt", header=T, sep="\t", row.names=1)

basal1 <- as.data.frame(basal$log2FoldChange, row.names=rownames(basal))
her21 <- as.data.frame(her2$log2FoldChange, row.names=rownames(her2))
luma1 <- as.data.frame(luma$log2FoldChange, row.names=rownames(luma))
lumb1 <- as.data.frame(lumb$log2FoldChange, row.names=rownames(lumb))
normal1 <- as.data.frame(normal$log2FoldChange, row.names=rownames(normal))

todos <- cbind(basal1, her21, luma1, lumb1, normal1)
sinormal <- cbind(basal1, her21, luma1, lumb1)















