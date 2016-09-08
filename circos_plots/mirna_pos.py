import csv

mirna_pos = {}

with open('mirna_pos') as f:
    r=csv.reader(f, delimiter=";")
    for l in r:
        # 'hsa-let-7a-1', 'MIMAT0004481', 'MIMAT0000062', 'chr9:94175957-94176036'
        (hsa, mimat1, mimat2, pos) = l

        pos = pos.replace("chr", "hs").replace(":"," ").replace("-"," ")
        mirna_pos[mimat1] = pos
        if mimat2:
            mirna_pos[mimat2] = pos


#from pprint import pprint
#pprint(mirna_pos)
