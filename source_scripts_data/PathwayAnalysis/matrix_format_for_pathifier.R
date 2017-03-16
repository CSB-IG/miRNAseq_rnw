# matrix for pathifier:
#----------------------------------------------------------------
# expression_rnaseq_tumour.txt: a matrix containing cases RNA-Seq expression data
# expression_rnaseq_control.txt: a matrix containing control RNA-Seq expression data
#-----------------------------------------------------------------

################################################
# Load expression data                         #
# Tab delimited non quoted .txt file matrices  #
# with the following format:                   #
#                                              #
# Name	Sample1	Sample2	Sample3	...            #
# Gene1	Value	Value	Value	...            #
# Gene2	Value	Value	Value	...	           #
# Gene3	Value	Value	Value	...            #
#                                              #
# The gene identifiers from the matrix should  #
# match the .gmt gene identifiers.             #
# Entrez gene ids are recommended              #
################################################

tumours <- read.table("expression_rnaseq_tumour.txt", header=T, row.names=1, sep = '\t')
controls <- read.table("expression_rnaseq_control.txt", header=T, row.names=1, sep = '\t')

# be sure that both matrices have the same genes in the same order
tumours <- tumours[intersect(rownames(tumours), rownames(controls)),]
controls <- controls[intersect(rownames(tumours), rownames(controls)),]

# assign code for controls and cases
code_cases <- rep("0", length(colnames(tumours)))
code_controls <- rep("1", length(colnames(controls)))

tumours <- rbind(code_case,tumours)
controls <- rbind(code_controls, controls)

# create matrix for pathifier in 
matrix <- cbind(tumours, controls)

write.table(data.frame("NORMALS"=rownames(matrix),matrix), "matrix_for_pathifier.txt", row.names=FALSE, sep="\t", quote=F)
