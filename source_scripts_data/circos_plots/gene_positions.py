import csv
from pprint import pprint

r = csv.reader(open('chr14/data/symbol2entrez.tsv'), delimiter=" ")
symbol = {}

for g in r:
    symbol[g[1]] = g[0]

r = csv.DictReader(open('chr14/data/mart_export.txt'), delimiter="\t")

positions = {}
for g in r:
    if g['chromosome'] =='X' \
       or g['chromosome'] == 'Y' \
       or g['chromosome'] in [str(a) for a in range(1,23)]:
        positions[symbol[g['entrezid']]] = "hs%s %s %s" % (g['chromosome'], g['start'], g['end'])
