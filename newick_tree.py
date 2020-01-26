#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:33:17 2019

@author: balthasar
"""


color = {
"ERR1100909" : "orange",
"ERR1100911" : "orange",
"ERR1100912" : "orange",
"ERR1100913" : "orange",
"ERR1100923" : "orange",
"ERR1100924" : "orange",
"ERR1100925" : "orange",
"ERR1100926" : "orange",
"ERR1100928" : "orange",
"ERR1100930" : "orange",
"ERR1100935" : "orange",
"ERR1100936" : "orange",
"ERR1100937" : "orange",
"ERR1100938" : "orange",
"ERR1100939" : "orange",
"ERR1100940" : "orange",
"ERR1100941" : "orange",
"ERR1100942" : "orange",
"ERR1100944" : "orange",
"ERR1100946" : "orange",
"ERR1100947" : "orange",
"ERR1100949" : "orange",
"ERR1100950" : "orange",
"ERR1100952" : "orange",
"ERR1100954" : "orange",
"ERR1100955" : "orange",
"ERR1100972" : "orange",
"ERR1100973" : "orange",
"ERR1100975" : "orange",
"ERR1100976" : "orange",
"0CEB535LM" : "orange",
"0CEB540LM" : "orange",
"0CEB550LM" : "orange",
"0CEB552LM" : "orange",
"0CEB553LM" : "orange",
"0CEB554LM" : "orange",
"0CEB559LM" : "orange",
"0CEB560LM" : "orange",
}


virulence = {
"ERR1100909" : 50.9,
"ERR1100911" : 71.3,
"ERR1100912" : 37.1,
"ERR1100913" : 61.5,
"ERR1100923" : 41.9,
"ERR1100924" : 58.1,
"ERR1100925" : 40.7,
"ERR1100926" : 31,
"ERR1100928" : 25,
"ERR1100930" : 50,
"ERR1100935" : 52.2,
"ERR1100936" : 65.5,
"ERR1100937" : 46.2,
"ERR1100938" : 65.5,
"ERR1100939" : 71.3,
"ERR1100940" : 45,
"ERR1100941" : 45.2,
"ERR1100942" : 61.5,
"ERR1100944" : 52.2,
"ERR1100946" : 29.3,
"ERR1100947" : 7,
"ERR1100949" : 68.2,
"ERR1100950" : 16,
"ERR1100952" : 18.2,
"ERR1100954" : 61.5,
"ERR1100955" : 40,
"ERR1100972" : 16,
"ERR1100973" : 16,
"ERR1100975" : 16,
"ERR1100976" : 16,
"0CEB535LM"  : 7,
"0CEB540LM"  : 7,
"0CEB550LM"  : 50.9,
"0CEB552LM"  : 50.9,
"0CEB553LM"  : 45,
"0CEB554LM"  : 45,
"0CEB559LM"  : 33.3,
"0CEB560LM"  : 33.3,
}

color_arla_vs_train = False

newv = {sn : ("blue" if v<25 else ("green" if v < 50 else ("orange" if v < 60 else "red")) ) for sn, v in virulence.items()}

Arla_results = [2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 3, 4, 2, 3, 4, 2, 2, 2, 2, 3, 3, 2, 2, 1, 3, 3,
 3, 2, 2, 2, 3, 2, 3, 3, 3, 3, 2, 2, 3]

colordict = {1:"blue", 2:"green",3:"orange",4:"red"}

arla_class_color = {str(i):colordict[c] for i, c in enumerate(Arla_results,1)}


if color_arla_vs_train:
    arla_class_color = {str(i):"blue" for i, c in enumerate(Arla_results,1)}
    newv = {sn : "red" for sn, v in virulence.items()}


color = {**arla_class_color, **newv }

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from Bio import Phylo


def label_func(string):
    string = str(string)
    if string[:3] == "DTU":
        return string[30:]
    else:
        return string

mpl.rcParams.update({'font.size': 5.5})

fname = "data/arla_accessory_binary_genes.fa.newick"


fig = plt.figure(figsize=(12, 7), frameon=False)
ax = fig.add_subplot(111)



patches = [mpatches.Patch(color=c, label='class ' + str(i)) for i, c in enumerate(["blue", "green", "orange", "red"],1)]
ax.legend(handles=patches)

tree = Phylo.read(fname, "newick")
tree.root_with_outgroup({'name': '24'})

Phylo.draw(tree, label_func=label_func, label_colors=color, axes=ax, show_confidence=False)



if color_arla_vs_train:
    patches = [mpatches.Patch(color=c, label=i) for i, c in zip(["Arla", "Training"], ["blue", "red"])]
    ax.legend(handles=patches)
    fig.savefig("plots/tree_arla_vs_train.pdf", bbox_inches = 'tight', pad_inches = 0)
else:
    patches = [mpatches.Patch(color=c, label='class ' + str(i)) for i, c in enumerate(["blue", "green", "orange", "red"],1)]
    ax.legend(handles=patches)
    fig.savefig("plots/tree_with_arla.pdf", bbox_inches = 'tight', pad_inches = 0)
    










