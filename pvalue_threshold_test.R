
enfermos <- data.matrix(read.table("~/TCGA_miRNA_BC/miRNAs_resultados/subredes/enfermos_p1_full_adjmtx.txt"))
colnames(enfermos) <- rownames(enfermos)
dim(enfermos)
isSymmetric(enfermos)

mir <- rownames(enfermos)[grepl("hsa*", rownames(enfermos))]
genes <- rownames(enfermos)[grepl("^(?!hsa*)", rownames(enfermos), perl=T)]

num_int_genes <- ((length(genes)**2)-length(genes))/2
num_int_mir <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 

(100/num_int_genes)*25334
100-((100/num_int_genes)*25334)

(100/num_int_mir)*14892
100-((100/num_int_mir)*14982)

# > (100/num_int_genes)*25334
# [1] 0.02211771
# > 100-((100/num_int_genes)*25334)                                                                                                                                                                                                                                                                                              
# [1] 99.97788
# > num_int_mir <- (length(mir)*length(genes))+(((length(mir)**2)-length(mir))/2) 
# > (100/num_int_mir)*14892
# [1] 0.1522526
# > 100-((100/num_int_mir)*14892)
# [1] 9.84775

# quantiles for 2533 miR and 1489 gene int

(100/num_int_mir)*2533
100-((100/num_int_mir)*2533)

(100/num_int_genes)*1489
100-((100/num_int_genes)*1489)

# > (100/num_int_mir)*2533
# [1] 0.02589684
# > 100-((100/num_int_mir)*2533)
# [1] 99.9741
# > 
# > (100/num_int_genes)*1489
# [1] 0.001299963
# > 100-((100/num_int_genes)*1489)
# [1] 99.9987

# quantiles for 253340 miR and 148920 gene int

(100/num_int_mir)*253340
100-((100/num_int_mir)*253340)

(100/num_int_genes)*148920
100-((100/num_int_genes)*1498920)

# > (100/num_int_mir)*253340
# [1] 2.590093
# > 100-((100/num_int_mir)*253340)
# [1] 99.77882
# > 
# > (100/num_int_genes)*148920
# [1] 0.1300138
# > 100-((100/num_int_genes)*148920)
# [1] 99.86999



