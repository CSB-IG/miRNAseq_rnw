# based on https://github.com/AngelCampos/Pathifier-Tool-Kit
# awk '{if(max<NF) {max=NF}} END {print max}' kegg_path_ch14.gmt
# awk '{ for( i=1; i<391; i++) {printf( $i "\t")}; printf($391 "\n")} ' kegg_path_ch14.gmt > prueba.gmt

# cut -f3- Project_wg_result1486230283_GSEA.gmt > 2

# printf 'NaN\n%.0s' {1..354} > 3

# paste -d"\t" pathways_names.txt NA pathways_genes.txt > wiki_path.gmt

################################################################################
# Running PATHIFIER (Drier et al., 2013)
# Author: Miguel Angel Garcia-Campos - Github: https://github.com/AngelCampos
################################################################################

# Install/load package ##############

library(pathifier)

# Load expression data for PATHIFIER
exp.matrix <- read.delim(file =file.choose(), as.is = T, row.names = 1)
# matriz_path_ensbl.txt

samples <- scan(file="sample_names_86_matched.txt", what="")
exp.matrix <- exp.matrix[,samples]

exp.matrix <- log2(exp.matrix+1)

# Load Genesets annotation 
gene_sets <- as.matrix(read.delim(file = file.choose(), header = F, sep = "\t", as.is = T))
                                  
# path.gmt

# filtra pathways con menos de 5 genes
gene_sets <- gene_sets[rowSums(is.na(gene_sets)) <= length(colnames(gene_sets))-5, ]

# filtra pathways con mas genes que samples
gene_sets <- gene_sets[rowSums(is.na(gene_sets)) >= length(colnames(gene_sets))-length(colnames(exp.matrix)), ]


#  Generate a list that contains genes in genesets
gs <- list()
for (i in 1:nrow(gene_sets)){
  a <- as.vector(gene_sets[i,3:ncol(gene_sets)])
  a <- na.omit(a)
  a <- a[a != ""]
  a <- matrix(a, ncol = 1)
  gs[[length(gs)+1]] <- a
  rm(a,i)
}

# Generate a list that contains the names of the genesets used
pathwaynames <- as.list(gene_sets[,1])

# Generate a list that contains the previos two lists: genesets and their names
PATHWAYS <- list(); PATHWAYS$gs <- gs; PATHWAYS$pathwaynames <- pathwaynames

# Prepare data and parameters ##################################################
# Extract information from binary phenotypes. 1 = Normal, 0 = Tumor
normals <- as.vector(as.logical(exp.matrix[1,]))
exp.matrix <- as.matrix(exp.matrix[-1, ])

# Calculate MIN_STD
N.exp.matrix <- exp.matrix[,as.logical(normals)]
rsd <- apply(N.exp.matrix, 1, sd)
min_std <- quantile(rsd, 0.25)

# Calculate MIN_EXP 
min_exp <- quantile(as.vector(exp.matrix), 0.1) # Percentile 10 of data

# Filter low value genes. At least 10% of samples with values over min_exp
# Set expression levels < MIN_EXP to MIN_EXP
over <- apply(exp.matrix, 1, function(x) x > min_exp)
G.over <- apply(over, 2, mean)
G.over <- names(G.over)[G.over > 0.1]
exp.matrix <- exp.matrix[G.over,]
exp.matrix[exp.matrix < min_exp] <- min_exp

# Set maximum 5000 genes with more variance
V <- names(sort(apply(exp.matrix, 1, var), decreasing = T))[1:5000]
V <- V[!is.na(V)]
exp.matrix <- exp.matrix[V,]
genes <- rownames(exp.matrix) # Checking genes
allgenes <- as.vector(rownames(exp.matrix))

# Generate a list that contains previous data: gene expression, normal status,
# and name of genes
DATASET <- list(); DATASET$allgenes <- allgenes; DATASET$normals <- normals
DATASET$data <- exp.matrix

# Run Pathifier ch14
PDS <- quantify_pathways_deregulation(DATASET$data, 
                                      DATASET$allgenes,
                                      PATHWAYS$gs,
                                      PATHWAYS$pathwaynames,
                                      DATASET$normals, 
                                      maximize_stability = T,
                                      attempts = 100,
                                      logfile="logfile_chr14.txt",
                                      min_std = min_std,
                                      min_exp = min_exp)

# Run Pathifier mir-200
PDS <- quantify_pathways_deregulation(DATASET$data, 
                                      DATASET$allgenes,
                                      PATHWAYS$gs,
                                      PATHWAYS$pathwaynames,
                                      DATASET$normals, 
                                      maximize_stability = T,
                                      attempts = 100,
                                      logfile="logfile_200.txt",
                                      min_std = min_std,
                                      min_exp = min_exp)
# Remove unnecesary data
rm(gene_sets, exp.matrix, allgenes, DATASET, PATHWAYS, rsd, V, over, G.over, 
   N.exp.matrix, gs, genes, min_exp, min_std, pathwaynames)

# Save image into working directory
save.image("PDS.RData")
