from mirna_pos import mirna_pos

from pprint import pprint

with open('mirnas.txt') as f:
    l = f.readlines()
    for mirna in l:
        (hsa,mimat) = mirna.split('.')
        mimat = mimat.strip()
                    
        print "%s;%s" % (mimat, mirna_pos.get(mimat, 'AGUAS con %s'%mimat))
