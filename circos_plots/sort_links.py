
import csv
from hsa_positions import position as hsa_pos
from mimat_positions import positions as mimat_pos

reader = csv.reader(open('/mnt/e/diana/TCGA_miRNA_BC/miRNAs_resultados/final/red.enfermos.MIR.86.txt'), delimiter="\t")

reader.next()
for l in reader:
    (s,t,w) = l
    if s.startswith('hsa') and t.startswith('hsa'):
        (s_hsa, s_mimat) = s.split('.') 
        (t_hsa, t_mimat) = t.split('.') 
        print s,t,hsa_pos.get(s_hsa, mimat_pos.get(s_mimat)),hsa_pos.get(t_hsa,mimat_pos.get(t_mimat))


# 2Tz2JhAe
