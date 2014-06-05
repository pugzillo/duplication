#!/usr/bin/python

####Python script that can take in a gene tree and detect paralogs and duplications; uses python module ete2
###linh c 6/5/2014

import sys

# Import ete2 module
from ete2 import PhyloTree

# Take in file input
file = sys.argv[1]
output = sys.argv[2]

# Load a tree structure from a newick file.
t = PhyloTree(file)

print t 

# Alternatively, you can scan the whole tree topology
events = t.get_descendant_evol_events()

# Open output file
fo = open(output, "wb+")

# Print its orthology and paralogy relationships
fo.write( 'Events detected from the root of the tree,' + file + '\n')
for ev in events:
    if ev.etype == "S":
        fo.write ('ORTHOLOGY RELATIONSHIP:' + ','.join(ev.in_seqs) + '<====>' + ','.join(ev.out_seqs) + '\n')
    elif ev.etype == "D":
        fo.write ('PARALOGY RELATIONSHIP:' + ','.join(ev.in_seqs) + '<====>' + ','.join(ev.out_seqs) + '\n')
fo.close()
