#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 18:41:38 2019

@author: balthasar
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


plt.close("all")

csv = np.loadtxt("data/all_genomes.csv", delimiter=",")
#arla_csv = np.loadtxt("supplemental_files/all_arla_genomes.csv", delimiter=",")

len_before_trimming = len(csv[0])
std_genes = np.std(csv, axis=0)
print(std_genes)

print(csv.shape)
std_sorted = np.sort(std_genes)
std_sumcount = np.arange(1,len(std_sorted)+1)
#
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.plot(std_sorted, std_sumcount, label='cumulative std count')
ax.axvline(15, color="red", alpha=0.3)
ax.axhline(2963, color="red", alpha=0.3)
ax.margins(x=0, y=0)
ax.set_xlabel("std")
ax.set_ylabel("cumulative count")
ax.text(14, 2000, "cutoff", rotation=90)
ax.text(10, 3000, "55.1%")
ax.set_xticks([0,10,15,20,30,40,50])
ax.set_yticks([0,1000,2000,2963,4000,5000])

fig.savefig("plots/std_count.pdf")

#arla_csv = arla_csv[:,np.std(csv, axis=0) > 15]
#csv = csv[:,np.std(csv, axis=0) > 15]



fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.imshow(csv, aspect="auto")
ax.set_xlabel("Gene")
ax.set_ylabel("Strain")
fig.colorbar(im, orientation="horizontal", fraction=0.1)
fig.savefig("plots/heatmap_blast.pdf", bbox_inches = 'tight', pad_inches = 0)
