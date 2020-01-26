#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 17:49:56 2020

@author: balthasar
"""

from ete3 import Tree, TreeStyle

strnr = [str(i) for i in range(51)]

t = Tree("data/arla_accessory_binary_genes.fa.newick")

t.set_outgroup("24")

ts = TreeStyle()
ts.show_leaf_name = True
ts.scale =  800

for n in t.traverse():
    if n.is_leaf():
        if n.name in strnr:
            n.img_style["fgcolor"] = "red"
        n.img_style["size"] = 10

t.render("plots/tree_ete.pdf", w=300, tree_style=ts)



t2 = Tree("data/accessory_binary_genes.fa.newick")

ts2 = TreeStyle()
ts2.show_leaf_name = False
ts2.mode = "c"
ts2.arc_start = -180 # 0 degrees = 3 o'clock
ts2.arc_span = 359


for n in t2.traverse():
    if n.is_leaf():
        if "DTU" in n.name:
            print(n.name)
            n.img_style["bgcolor"] = "red"
        if "CEB" in n.name or "ERR" in n.name:
            print(n.name)
            n.img_style["bgcolor"] = "blue"

t2.render("plots/tree_ete2.png", w=1200, tree_style=ts2)

