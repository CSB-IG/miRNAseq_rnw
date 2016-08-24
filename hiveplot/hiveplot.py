from pyveplot import *
import networkx as nx
import csv
import argparse


parser = argparse.ArgumentParser(description='hiveplot of mirna-gene interactions')
parser.add_argument('edgelist', type=argparse.FileType('r'), help="network in edgelist format")
args   = parser.parse_args()



g = nx.Graph()

for edge in csv.DictReader(args.edgelist, delimiter="\t"):
    g.add_edge( edge['from'], edge['to'], w=edge['weight'] )


by_degree = {}
for n in g.nodes():
    d = g.degree(n)
    if d in by_degree:
        by_degree[d].append(n)
    else:
        by_degree[d] = [n, ]

    
# nx.write_gpickle(g, "test.gpickle")

from pprint import pprint
pprint(by_degree)
