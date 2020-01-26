#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 16:03:27 2019

@author: balthasar
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


#Core genes      (99% <= strains <= 100%)        881
#Soft core genes (95% <= strains < 99%)  658
#Shell genes     (15% <= strains < 95%)  2237
#Cloud genes     (0% <= strains < 15%)   39509
#Total genes     (0% <= strains <= 100%) 43285

plt.close("all")



ind = np.arange(3)
width = 0.35


for run in range(1,4):

    genes = {
    "core genes      (.99 -   1)" : np.array([881,   1200, 2191]),
    "soft core genes (.95 - .99)" : np.array([658,   1022, 39  ]),
    "shell genes     (.15 - .95)" : np.array([2237,  1352, 1297]),
    "cloud genes     (  0 - .15)" : np.array([39509, 3869, 1847]),
    }

    labels = np.array(('NCBI strains', 'Arla strains', 'Training strains'))
    bottom = np.array((0,0,0))
    
    fig = plt.figure()
    ax2, ax = fig.subplots(2, 1,sharex=True)
    fig.subplots_adjust(hspace = 0.05)
    
    for gtype, count in genes.items():
        count[run:] = 0
        ax.bar(ind, count, width,bottom = bottom, label=gtype)
        ax2.bar(ind, count, width,bottom = bottom, label=gtype)
        for i in ind:
            if count[i]:
                ax.text(i, bottom[i] + (count[i] if count[i] < 4000 else 4000)/2, str(count[i]),
                        ha='center', va='center', color = 'white',
                        path_effects = [path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        bottom += count
    
    
    
    ax.set_ylim([0, 8300])
    ax.set_ylabel('Gene count')
    ax.set_xticks(ind)
    labels[run:] = ""
    ax.set_xticklabels(labels)
    
    #ax.text()
    
    #plt.yticks(np.arange(0, 81, 10))
    ax.spines['top'].set_visible(False)
    
    
    ax2.spines['bottom'].set_visible(False)
    ax2.tick_params(bottom=False)
    ax2.set_ylim([39700, 45000])
    ax2.legend(prop={'family': 'DejaVu Sans Mono'})
    d = .015
    
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
    ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    
    kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
    ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    fig.savefig(f"plots/barchart_{str(run)}.pdf")