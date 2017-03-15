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


h = Hiveplot( args.plot.name )
args.plot.close()
h.dwg.width=15000
h.dwg.height=15000

centre = (750, 500)
centre1 = (750-15, 500+20)


# configure mirna axes
m1 = Axis( coords(20, -30 -115, centre), coords(500, -30 -115, centre), stroke="#CBFF65", stroke_width=10)
m2 = Axis( coords(20, -30, centre), coords(500, -30, centre), stroke="#CBFF65", stroke_width=10)
pos   = 0.0
delta = 1.0 / len(mirnas)
for n in mirnas:
    node0 = Node(n)    
    node1 = Node(n)
    m1.add_node( node0, pos ) 
    m2.add_node( node1, pos )
    pos += delta


g1 = Axis(coords(20, 90, centre), coords(500, 90, centre), stroke="#00C598", stroke_width=10)
g2 = Axis(coords(20, -30 -115, centre1), coords(500, -30-115, centre1), stroke="#00C598", stroke_width=10)

pos   = 0.0
delta = 1.0 / len(genes)
for n in genes:
    node0 = Node(n)
    node1 = Node(n) 
    g1.add_node( node0, pos ) 
    g2.add_node( node1, pos )
    pos += delta

h.axes = [m1, m2, g1, g2]


bar = progressbar.ProgressBar()

for e in bar(g.edges()):

    if e[0] in mirnas and e[1] in mirnas:     # mirnas to mirnas
        h.connect(m1, e[0], 45,
                  m2, e[1], -45,
                  stroke_width   = g.get_edge_data(*e)['w'] * 4,
                  stroke_opacity = 0.1,
                  stroke         = 'grey')

    if e[0] in mirnas and e[1] in genes: # mirnas to genes
        h.connect(m2, e[0], 45,
                  g1, e[1], -45,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=1/100.0,
                  stroke="grey")

    if e[1] in mirnas and e[0] in genes: # mirnas to genes
        h.connect(m2, e[1], 45,
                  g1, e[0], -45,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=1/100.0,
                  stroke="grey")

        

    if e[0] in genes and e[1] in genes: # genes to genes
        h.connect(g1, e[0], 45,
                  g2, e[1], -45,
                  stroke_width   = g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity = 0.06,
                  stroke='grey')

print "saving"        
h.save()
