import argparse
import csv
from hsa_positions import position as hsa_pos
from mimat_positions import positions as mimat_pos

parser = argparse.ArgumentParser(description='Create a links file for circos plot from edgelist.')
parser.add_argument('edgelist', type=argparse.FileType('r'), help="path to edgelist" )
args = parser.parse_args()

reader = csv.reader(args.edgelist, delimiter=";")

reader.next()
for l in reader:
    (s,t,w) = l
    if s.startswith('hsa') and t.startswith('hsa'):
        (s_hsa, s_mimat) = s.split('.') 
        (t_hsa, t_mimat) = t.split('.') 
        print hsa_pos.get(s_hsa, mimat_pos.get(s_mimat)), hsa_pos.get(t_hsa,mimat_pos.get(t_mimat))


# 2Tz2JhAe
