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
    

genes       = []
mirnas      = []
# classify nodes
for n in degree_ordered:
    if not n.startswith('hsa'):
        genes.append(n)
    else:
        mirnas.append(n)



h = Hiveplot( args.plot.name )
args.plot.close()

h.dwg.width=15000
h.dwg.height=15000



# configure mirna axes
m1 = Axis( (20,750), (20, 20), stroke="maroon", stroke_width=4)
m2 = Axis( (350,750), (350, 20), stroke="maroon", stroke_width=4)
pos   = 0.0
delta = 1.0 / len(mirnas)
for n in mirnas:
    node0 = Node(n)    
    node1 = Node(n)
    m1.add_node( node0, pos ) 
    m2.add_node( node1, pos )
    pos += delta


# genes axes
g1 = Axis( (680, 1000), (680, 20), stroke="dodgerblue", stroke_width=4)
g2 = Axis( (1010, 1000), (1010, 20), stroke="dodgerblue", stroke_width=4)
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
        h.connect(m1, e[0], 5,
                  m2, e[1], -5,
                  stroke_width   = g.get_edge_data(*e)['w'] * 2,
                  stroke_opacity = 0.5,
                  stroke         = 'burlywood')


    if e[0] in mirnas and e[1] in genes:
        h.connect(m2, e[0], 5,
                  g1, e[1], -5,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=0.09,
                  stroke='coral')

    if e[1] in mirnas and e[0] in genes:
        h.connect(m2, e[1], 5,
                  g1, e[0], -5,
                  stroke_width=g.get_edge_data(*e)['w'] * 5,
                  stroke_opacity=0.09,
                  stroke='coral')

    if e[0] in genes and e[1] in genes:
        h.connect(g1, e[0], 5,
                  g2, e[1], -5,
                  stroke_width=g.get_edge_data(*e)['w'] * 6,
                  stroke_opacity=0.04,
                  stroke='darksalmon')

print "saving"        
h.save()
