import csv
import argparse
import networkx as nx

parser = argparse.ArgumentParser(description='Get first neighbors of region.')
parser.add_argument('--edgelist', type=argparse.FileType('r'), help="path to edgelist" )
parser.add_argument('--chrom', type=str, help="chromosome" )
parser.add_argument('--start', type=int, help="start position" )
parser.add_argument('--end', type=int, help="end position" )
parser.add_argument('--output', type=argparse.FileType('w'), help="filtered edgelist" )
args = parser.parse_args()


G = nx.Graph()

# create network from edgelist
reader = csv.reader(args.edgelist, delimiter=" ")

for l in reader:
    s = tuple(l[0:3])
    t = tuple(l[3:6])
    w = l[6].split('=')[1]
    G.add_edge( s, t, w = w )


H = nx.Graph()
for n in G.nodes():
    if n[0] == args.chrom and int(n[1])>=args.start and int(n[1])<=args.end:
        for m in G.neighbors(n):
            H.add_edge(n, m, w=G[n][m]['w'])

            
o=csv.writer(args.output, delimiter=" ")
for e in H.edges():
    row = list(e[0]) + list(e[1]) + ["w=%s" % H.get_edge_data(*e)['w'],]
    o.writerow(row)


