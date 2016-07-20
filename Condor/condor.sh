

# 0) source ~/enviromets/csbig/bin/activate
# 1) modificar expfile, probes, run_id, y outdir
#probes no debe de tener ningun header

python /home/diana/breast_cancer_networks/parallel_aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/enfermos.86.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/GENEID_filtro.txt \
    --run_id RNAseq_86.enfermos.rsem \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/datos_nivel3_norm_TCGA/RNAseq_86.enfermos.rsem \
    --p 1

# 2)  Poner en el encabezado de .condor abajo de log
# requirements = Machine == "notron.inmegen.gob.mx"

# 3) Dentro del outdir darle en linea de comandos:
# condor_submit RNAseq_86.enfermos.rsem.condor

# 4) condor_status

sanos <- "sanos.86.txt"
probes sanos <- "GENEID.86.txt"

Basal
python /home/diana/rnw/parallel_aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_Basal_Expmtx_Aracne_RNAseq.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_RNAseq_752 \
    --run_id RNAseq_752_casos_red_Basal \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/RNAseq_752_casos_red_Basal \
    --p 1
    
python /home/diana/parallel-aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_Basal_Expmtx_Aracne.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_miRNAs_752.txt \
    --run_id miRNAseq_752_casos_red_Basal \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/miRNAseq_752_casos_red_Basal \
    --p 1

LumA
python /home/diana/rnw/parallel_aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_LumA_Expmtx_Aracne_RNAseq.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_RNAseq_752 \
    --run_id RNAseq_752_casos_red_LumA \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/RNAseq_752_casos_red_LumA \
    --p 1
    
python /home/diana/parallel-aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_LumA_Expmtx_Aracne.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_miRNAs_752.txt \
    --run_id miRNAseq_752_casos_red_LumA \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/miRNAseq_752_casos_red_LumA \
    --p 1

LumB
python /home/diana/rnw/parallel_aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_LumB_Expmtx_Aracne_RNAseq.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_RNAseq_752 \
    --run_id RNAseq_752_casos_red_LumB \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/RNAseq_752_casos_red_LumB \
    --p 1
    
python /home/diana/parallel-aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_LumB_Expmtx_Aracne.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_miRNAs_752.txt \
    --run_id miRNAseq_752_casos_red_LumB \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/miRNAseq_752_casos_red_LumB \
    --p 1

Her2
python /home/diana/rnw/parallel_aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_Her2_Expmtx_Aracne_RNAseq.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_RNAseq_752 \
    --run_id RNAseq_752_casos_red_Her2 \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/RNAseq_752_casos_red_Her2 \
    --p 1
    
python /home/diana/parallel-aracne/genera_condor.py \
    --path_to_aracne2 /home/diana/ARACNE/aracne2 \
    --expfile /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/Subtipo_Her2_Expmtx_Aracne.txt \
    --probes /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/probes_miRNAs_752.txt \
    --run_id miRNAseq_752_casos_red_Her2 \
    --outdir /mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/miRNAseq_752_casos_red_Her2 \
    --p 1
