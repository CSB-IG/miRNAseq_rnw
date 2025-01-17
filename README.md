# miRNAseq_rnw
miRNAseq regulatory networks


1. expression matrix mature miRNA (miRBase_mysql.txt):
> dim(x)
[1] 2588  752
2. formating matrix (format_miR_matrix.R):
matrix format: rownames: GENE/miR ID, colnames: patient code: XXXX.01(case)/11(control)
> miRNAs_752_maduros.txt
3. filter raw counts: filter targets with a minimum of  5  counts  in  at  least  25%  of  the  samples(752) = 188
low_counts_filter_miRNA.R
 > dim(miRNA)
 [1] 715 752
miRNAs_752_maduros_filtro.txt
4. density plot
> density_plot.miR_RNAseq.R
5. Normalization
miRNA_normalization.R
"Adjusting  the data by cpm or total count scaling introduces more variability to
the  data,  whereas  all  other  methods  resulted  in  more  similar
distribution  across  all  samples."
" The distributions are more similar and centered
at zero when the data are normalized using UQ, TMM, DESeq,
cyclic loess and quantile normalizations"
"UQ  and  TMM  decreased  the  variance  across  all miRNAs compared with the raw data"
R statistical environment (v3.1.2)
Bioconductor  (3.0)
Cyclic loess and quantile normalization were performed using  the normalizeBetweenArrays() function  in  the  limma package (v3.22.4)
TMM and UQ normalizations were performed using calcNormFactors() in the edgeR (v3.8.5) package
the  function estimateSizeFactors() in  the  DESeq (v1.18.0) package was used to normalize the count data using size
factors 
(http://bib.oxfordjournals.org/content/early/2015/04/17/bib.bbv019.full.pdf+html)
> miRNA_752_normTMM.txt
> TMM1.miRNA.753.pdf
> raw.miRNA.753.pdf
6. add RNAseq expression data to miRNAseq for mi miRNA-miRNA miRNA-mRNA interactions
expression_matrix_miRNAs.R
7. add RNAseq network mRNA-mRNA interaction for whole network
> red.RNAseq.final.752.txt = 99.987% dpi 0.1 11202 interactions (with no header to cat under miRNA network)
> red.MIRseq.final.752.txt = 99.980% dpi 0.1 
cat red.MIRseq.final.752.txt red.RNAseq.final.752.txt > red.final.752.txt
>>> red.final.752.txt

for 86 samples matrices
enfermos <- read.table("miRNAs_enfermos_maduros.txt", sep="\t", header=T, row.names=1)
sanos <- read.table("miRNAs_sanos_maduros.txt", sep="\t", header=T, row.names=1)
1. expression matrix mature miRNA (miRBase_mysql.txt):
> dim(enfermos)
[1] 2588  752
> dim(sanos)
[1] 2588  752
2. formating matrix (format_miR_matrix.R):
matrix format: rownames: GENE/miR ID, colnames: patient code: XXXX.01(case)/11(control)
> miRNAs_enfermos_maduros.txt
> miRNAs_sanos_maduros.txt
3. filter raw counts: filter targets with a minimum of  5  counts  in  at  least  25%  of  the  samples(86) = 22
low_counts_filter_miRNA.R
> dim(enfermos)                                                                                                                       
[1] 690  86
> dim(sanos)
[1] 658  86
> miRNAs_enfermos_maduros_filtro.txt
> miRNAs_sanos_maduros_filtro.txt
set two groups to intersection (633)
miRNAs_752_maduros_filtro.txt
4. density plot
> density_plot.miR_RNAseq.R
5. Normalization
miRNA_normalization.R
"Adjusting  the data by cpm or total count scaling introduces more variability to
the  data,  whereas  all  other  methods  resulted  in  more  similar
distribution  across  all  samples."
" The distributions are more similar and centered
at zero when the data are normalized using UQ, TMM, DESeq,
cyclic loess and quantile normalizations"
"UQ  and  TMM  decreased  the  variance  across  all miRNAs compared with the raw data"
R statistical environment (v3.1.2)
Bioconductor  (3.0)
Cyclic loess and quantile normalization were performed using  the normalizeBetweenArrays() function  in  the  limma package (v3.22.4)
TMM and UQ normalizations were performed using calcNormFactors() in the edgeR (v3.8.5) package
the  function estimateSizeFactors() in  the  DESeq (v1.18.0) package was used to normalize the count data using size
factors 
(http://bib.oxfordjournals.org/content/early/2015/04/17/bib.bbv019.full.pdf+html)
> miRNAs_86_maduros_norm.txt
> miRNAs_86_maduros_norm_enfermos.txt
> miRNAs_86_maduros_norm_sanos.txt
> TMM.miRNA.86.pdf
> raw.miRNA.86.pdf
6. add RNAseq expression data to miRNAseq
expression_matrix_miRNAs.R
7. add RNAseq network mRNA-mRNA interaction for whole network
> enfermos.99.987.dpi.0.1.sif.txt -> red.enfermos.RNAseq.86.txt = 99.987% dpi 0.1 14280 interactions (with no header to cat under miRNA network)
> sanos.99.987.dpi.0.1.sif.txt > red.sanos.RNAseq.86.txt = 99.987% dpi 0.1 14890 interactions
> sanos.mir.99.980.dpi.sif.txt > red.sanos.MIR.86.txt = 99.980% dpi 0.1 24401 interactions
cat red.sanos.MIR.86.txt red.sanos.RNAseq.86.txt > red.final.sanos.txt
> enfermos.mir.99.980.dpi.sif.txt > red.enfermos.MIR.86.txt = 99.980% dpi 0.1 23582 interactions
cat red.enfermos.MIR.86.txt red.enfermos.RNAseq.86.txt > red.final.enfermos.txt

probes
* miR752 -> probes_miRNA.txt
yes 0 | head -n 715 > 1.txt
paste probes_miRNA.txt 1.txt > clave_MIR_752.txt
* RNAseq752 -> probes_RNAseq_752
yes 1 | head -n 15972 > 1.txt
paste probes_RNAseq_752 1.txt > clave_RNA_752.txt
* miR86 -> probes_86.txt
yes 0 | head -n 633 > 1.txt
paste probes_86.txt 1.txt > clave_MIR_86.txt
*RNAseq86 -> probes_RNAseq_86
yes 1 | head -n 15136 > 1.txt
paste probes_RNAseq_86 1.txt > clave_RNA_86.txt

cat clave_MIR_86.txt RNAseq_86_genes.txt RNAseq_86_TF.txt > clave.final.86.txt
cat clave_MIR_752.txt RNAseq_752_genes.txt RNAseq_752_TF.txt > clave.final.752.txt

enfermos <- read.table(file = "enfermos.86.adj.txt", header = T, sep = '\t')
p <- enfermos
p[p < 0.16482207] <- 0
sum(enfermos != 0)/2
sum(p != 0)/2

