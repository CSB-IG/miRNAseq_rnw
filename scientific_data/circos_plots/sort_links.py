import argparse
import csv

from gene_positions import positions as gene_pos
from mirna_pos import mirna_pos

parser = argparse.ArgumentParser(description='Create a links file for circos plot from edgelist.')
parser.add_argument('edgelist', type=argparse.FileType('r'), help="path to edgelist" )
args = parser.parse_args()

reader = csv.DictReader(args.edgelist, delimiter="\t")

def node2pos(node):
    if node.startswith('hsa-'):
        (hsa, mimat) = node.split('.')
        return mirna_pos.get(mimat, "AGUAS")
    else:
        return gene_pos.get(node)

for l in reader:
    if node2pos(l['from']) and node2pos(l['to']):
        print node2pos(l['from']), node2pos(l['to']), "w=%s" % l['weight']

