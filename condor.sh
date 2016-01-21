

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