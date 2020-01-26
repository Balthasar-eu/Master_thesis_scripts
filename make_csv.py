#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 08:32:33 2019

@author: balthasar
"""

# This script will make the table for blastp (table for heatmap)
# input file = list of query name , best hit from blastp results
# output file = table of present and absent of query (table for heatmap)
#### EXAMPLE ####
# python make_table_from_blastp_percentIdentity_16June15.py list_query_name  blast_result output

import sys

# open list of query name
filepath1 = sys.argv[1]
with open(filepath1,'r') as f:
    query = [line.strip() for line in f]

# open blast result file
filepath2 = sys.argv[2]
with open(filepath2, 'r') as f:
    blastout = [line.split('\t') for line in f]

hits = {}
for line in blastout:
    if line[0] in hits:
        if line[3] > hits[line[0]][1]:
             hits[line[0]] = [line[2],line[3]]
    else:
        hits[line[0]] = [line[2],line[3]]

filepath3 = sys.argv[3]
with open(filepath3,'w') as f:
    f.write(str(filepath2))
    for q in query:
        f.write(','+str(hits[q][0])) if q in hits else f.write(',0')
    f.write('\n')
