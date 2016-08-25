from pyveplot import *
import networkx as nx
import csv
import argparse
import math


def coords( radius, angle, origin=(0,0)):
    """ Returns a tuple of tuples of coordinates given a radius, angle and origin tuple """
    return (origin[0] + round((radius * math.cos(math.radians(angle))), 2),
            origin[1] + round((radius * math.sin(math.radians(angle))), 2))


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
for d in sorted(by_degree.keys(), reverse=False):
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
        genes.append(n)        
        if only_talks_to_mrna(n):
            mirna_genes.append(n)
    else:
        mirnas.append(n)


print len(mirnas),len(genes),len(mirna_genes)


h = Hiveplot( args.plot.name )
args.plot.close()
h.dwg.size=(5000,5000)

centre = (750, 250)


# configure mirna axes
m1 = Axis( centre, coords(250, -120, centre), stroke="red", stroke_width=4)
m2 = Axis( centre, coords(250, -60, centre), stroke="red", stroke_width=4)
pos   = 0.0
delta = 1.0 / len(mirnas)
for n in mirnas:
    node0 = Node(n)    
    node1 = Node(n)
    m1.add_node( node0, pos ) 
    m2.add_node( node1, pos )
    pos += delta



ga = Axis(centre, coords(750, 30, centre), stroke="orange", stroke_width=4)
mg = Axis(centre, coords(750, 150, centre), stroke="red", stroke_width=4)


pos   = 0.0
delta = 1.0 / len(genes)
for n in genes:
    node0 = Node(n)
    node1 = Node(n) 
    ga.add_node( node0, pos ) 
    mg.add_node( node1, pos )
    pos += delta

h.axes = [m1, m2, ga, mg]

    
for e in g.edges():
    # mirnas to mirnas
    if e[0] in mirnas and e[1] in mirnas:
        h.connect(m1, e[0], 25,
                  m2, e[1], -25,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='purple')


    if e[0] in mirnas and e[1] in mirna_genes: # mirnas to mirna-genes
        h.connect(m1, e[0], -25,
                  mg, e[1], 25,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='blue')
    elif e[0] in mirnas and e[1] in genes: # mirnas to genes
        h.connect(m2, e[0], 25,
                  ga, e[1], -25,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='green')


    if e[1] in mirnas and e[0] in mirna_genes: # mirna-genes to mirnas
        h.connect(m1, e[1], -25,
                  mg, e[0], 25,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='blue')
    elif e[1] in mirnas and e[0] in mirna_genes: # genes to mirnas
        h.connect(m2, e[1], 25,
                  ga, e[0], -25,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='green')


    if e[0] in genes and e[1] in genes or e[1] in genes and e[0] in genes: # gene to gene
        h.connect(ga, e[0], 25,
                  mg, e[1], -25,
                  stroke_width='0.34', stroke_opacity='0.4', stroke='yellow')

        
h.save()

# nx.write_gpickle(g, "test.gpickle")

#pprint(mirna_genes)
