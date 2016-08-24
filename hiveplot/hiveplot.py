from pyveplot import *
import networkx as nx
import csv
import argparse

from pprint import pprint

parser = argparse.ArgumentParser(description='hiveplot of mirna-gene interactions')
parser.add_argument('edgelist', type=argparse.FileType('r'), help="network in edgelist format")
parser.add_argument('plot', type=argparse.FileType('w'), help="plot file")
args   = parser.parse_args()




# load edges into graph
g = nx.Graph()
for edge in csv.DictReader(args.edgelist, delimiter="\t"):
    g.add_edge( edge['from'], edge['to'], w=edge['weight'] )


# sort nodes by degree
by_degree = {}
for n in g.nodes():
    d = g.degree(n)
    if d in by_degree:
        by_degree[d].append(n)
    else:
        by_degree[d] = [n, ]

degree_ordered = []
for d in sorted(by_degree.keys(), reverse=True):
    degree_ordered += by_degree[d]



# if a gene doesn't have other genes as its first neighbors
def only_talks_to_mrna( gene ):
    status = True
    for n in g.neighbors(gene):
        if not n.startswith('hsa'):
            return False
    return status
    


mirna_genes = []
genes       = []
mirnas      = []
# classify nodes
for n in degree_ordered:
    if not n.startswith('hsa'):
        if only_talks_to_mrna(n):
            mirna_genes.append(n)
        else:
            genes.append(n)
    else:
        mirnas.append(n)


print len(mirnas),len(genes),len(mirna_genes)


h = Hiveplot( args.plot.name )
args.plot.close()

a = Axis( (200,200), (200,100), stroke="grey" )
b = Axis( (250,200), (250,100), stroke="blue" )

h.axes = [a, b]

pos   = 0.0
delta = 1.0 / len(mirnas)
for n in mirnas:
    node0 = Node(n)    
    node1 = Node(n)
    a.add_node( node0, pos ) 
    b.add_node( node1, pos )
    pos += delta

for e in g.edges():
    if e[0] in mirnas and e[1] in mirnas:
        h.connect(a, e[0], 5,
                  b, e[1], -5,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='purple')
    # if e[0] in axis_mrna1.nodes and e[1] in axis_mrna0.nodes:
    #     h.connect(axis_mrna1, e[1], 20,
    #               axis_mrna0, e[0], -20,
    #               stroke_width='0.34', stroke_opacity='0.4', stroke='purple')
        
h.save()

# nx.write_gpickle(g, "test.gpickle")

#pprint(mirna_genes)
