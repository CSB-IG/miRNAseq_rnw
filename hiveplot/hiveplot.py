from pyveplot import *
import networkx as nx
import csv
import argparse
import math
import progressbar

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
    g.add_edge( edge['from'], edge['to'], w=float(edge['weight']) )


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
h.dwg.width=15000
h.dwg.height=15000

centre = (750, 500)


# configure mirna axes
m1 = Axis( coords(70, -100, centre), coords(500, -105, centre), stroke="maroon", stroke_width=4)
m2 = Axis( coords(70, -80, centre), coords(500, -75, centre), stroke="maroon", stroke_width=4)
pos   = 0.0
delta = 1.0 / len(mirnas)
for n in mirnas:
    node0 = Node(n)    
    node1 = Node(n)
    m1.add_node( node0, pos ) 
    m2.add_node( node1, pos )
    pos += delta



ga = Axis(coords(50, 45, centre), coords(750, 45, centre), stroke="dodgerblue", stroke_width=4)
mg = Axis(coords(50, 135, centre), coords(750, 135, centre), stroke="dodgerblue", stroke_width=4)


pos   = 0.0
delta = 1.0 / len(genes)
for n in genes:
    node0 = Node(n)
    node1 = Node(n) 
    ga.add_node( node0, pos ) 
    mg.add_node( node1, pos )
    pos += delta

h.axes = [m1, m2, ga, mg]


bar = progressbar.ProgressBar()

for e in bar(g.edges()):

    if e[0] in mirnas and e[1] in mirnas:     # mirnas to mirnas
        h.connect(m1, e[0], 5,
                  m2, e[1], -5,
                  stroke_width   = g.get_edge_data(*e)['w'] * 2,
                  stroke_opacity = 0.5,
                  stroke         = 'purple')


    if e[0] in mirnas and e[1] in mirna_genes: # mirnas to mirna-genes
        h.connect(m1, e[0], -45,
                  mg, e[1], 45,
                  stroke_width   = g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity = 0.09,
                  stroke         = 'slategray')
    elif e[0] in mirnas and e[1] in genes: # mirnas to genes
        h.connect(m2, e[0], 45,
                  ga, e[1], -45,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=0.09,
                  stroke='gainsboro')


    if e[1] in mirnas and e[0] in mirna_genes: # mirna-genes to mirnas
        h.connect(m1, e[1], -45,
                  mg, e[0], 45,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=0.09,
                  stroke='slategray')
    elif e[1] in mirnas and e[0] in mirna_genes: # genes to mirnas
        h.connect(m2, e[1], 45,
                  ga, e[0], -45,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=0.09,
                  stroke='gainsboro')


    if e[0] in genes and e[1] in genes or e[1] in genes and e[0] in genes: # gene to gene
        h.connect(ga, e[0], 33,
                  mg, e[1], -33,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=0.05,
                  stroke='tan')

print "saving"        
h.save()
