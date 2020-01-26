#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 19:29:57 2019

@author: balthasar
"""

import numpy as np
from sklearn.decomposition import PCA
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats.stats import pearsonr

from sklearn.svm import SVC
import sklearn

plt.close("all")

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
"L_monocytogenes_10CEB535LM" : 7,
"L_monocytogenes_10CEB540LM" : 7,
"L_monocytogenes_10CEB550LM" : 50.9,
"L_monocytogenes_10CEB552LM" : 50.9,
"L_monocytogenes_10CEB553LM" : 45,
"L_monocytogenes_10CEB554LM" : 45,
"L_monocytogenes_10CEB559LM" : 33.3,
"L_monocytogenes_10CEB560LM" : 33.3,
}


clinical_freq = np.array([50.9,71.3,37.1,61.5,41.9,58.1,40.7,31,25,50,52.2,65.5,
                          46.2,65.5,71.3,45,45.2,61.5,52.2,29.3,7,68.2,16,18.2,
                          61.5,40,16,16,16,16,7,7,50.9,50.9,45,45,33.3,33.3])

newv = [(1 if v<25 else (2 if v < 50 else (3 if v < 60 else 4)) ) for v in clinical_freq]

csv = np.loadtxt("data/all_genomes.csv", delimiter=",")

path = "kma_output/"

with open("data/query_list.txt") as f:
    querylist = [l[:-1] for l in f]
 
dirlist = os.listdir(path)

with open("data/strains.txt", "r") as f:
    strains = [l[:-1] for l in f]

templatelist = []

for d in strains:
    with open(path + d + ".res") as f:
        templatelist.append({line.split("\t")[0].split(" ")[0]:float(line.split("\t")[4]) for i, line in enumerate(f) if i})


final = np.array([[tdict[entry] if entry in tdict else 0 for tdict in templatelist] for entry in querylist]).T

np.savetxt("data/kma_genomes.csv",final,fmt="%.2f", delimiter=",")

p_kma_blast = []

for k,c in zip(final,csv):
    p_kma_blast.append(pearsonr(k,c)[0])

pearson_all = []

for i in range(len(final)-1,-1,-1):
    for j in range(i):
        pearson_all.append(pearsonr(final[i],final[j])[0])

print("KMA to blast:", np.mean(p_kma_blast))
print("All:", np.mean(pearson_all))

final = final[:,np.std(csv, axis=0) > 15]
csv = csv[:,np.std(csv, axis=0) > 15]

p_kma_blast = []

for k,c in zip(final,csv):
    p_kma_blast.append(pearsonr(k,c)[0])

pearson_all = []

for i in range(len(final)-1,-1,-1):
    for j in range(i):
        pearson_all.append(pearsonr(final[i],final[j])[0])

print("KMA to blast:", np.mean(p_kma_blast))
print("All:", np.mean(pearson_all))

pca = PCA(n_components=3)
pca.fit(np.vstack((csv,final)))

csvtrans = pca.transform(csv)
finaltrans = pca.transform(final)
alltrans = pca.transform(np.vstack((csv,final)))
###### plot 

fig = plt.figure(figsize=(8, 8), frameon=False)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')

fig.patch.set_visible(False)
ax.patch.set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])


ax.scatter(csvtrans[:,0],csvtrans[:,1],csvtrans[:,2], label="blast")
ax.scatter(finaltrans[:,0],finaltrans[:,1],finaltrans[:,2], label="kma")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels)
#ax.scatter(alltrans[:,0],alltrans[:,1],alltrans[:,2], c=2*newv)

for x,y,z,n in zip(csvtrans[:,0],csvtrans[:,1],csvtrans[:,2],range(38)):
    ax.text(x,y,z,n)

for x,y,z,n in zip(finaltrans[:,0],finaltrans[:,1],finaltrans[:,2],range(38)):
    ax.text(x,y,z,n)

fig.savefig('plots/compare_kma.pdf')



###############################################################################
############################ Compare ML ###############################
###############################################################################

csv = np.loadtxt("data/all_genomes.csv", delimiter=",")

std_genes = np.std(csv, axis=0)

std_sorted = np.sort(std_genes)
std_sumcount = np.arange(1,len(std_sorted)+1)

csv = csv[:,np.std(csv, axis=0) > 15]

newv = {sn: (1 if v<25 else (2 if v < 50 else (3 if v < 60 else 4)) ) for sn,v in virulence.items()}

count = {1:0,2:0,3:0,4:0}
for k in newv:
    count[newv[k]] += 1

with open("data/keyfile.txt","r") as f:
    keys = [line for line in f]

Y = np.array([newv[k[:-10]] for k in keys], dtype=int)

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=1)
kma_train, kma_test, kma_y_train, kma_y_test = sklearn.model_selection.train_test_split(final, Y, random_state=1)

clf = SVC(kernel="linear", C=1)
clf.fit(X_train, y_train)

print(clf.score(X_test,y_test))
print(clf.score(kma_test,y_test))

scores = []

for i in range(100):
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=i)
    kma_train, kma_test, kma_y_train, kma_y_test = sklearn.model_selection.train_test_split(final, Y, random_state=i)
    clf = SVC(kernel="linear", C=1)
    clf.fit(X_train, y_train)
    scores.append(clf.score(kma_test,y_test))


scores = np.array(scores)
values, counts = np.unique(scores, return_counts=True)

fig = plt.figure(frameon=False)
ax = fig.add_subplot(111)
ax.bar(values, counts, 0.05)
ax.set_title(f"Mean: {np.mean(scores):.2f} σ: {np.std(scores):.2f}")

ax.set_xlabel("score")
ax.set_ylabel("count")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels)
fig.savefig("plots/kma_scores.pdf")


###############################################################################
# KMA ONLY

scores = []

for i in range(100):
    kma_train, kma_test, kma_y_train, kma_y_test = sklearn.model_selection.train_test_split(final, Y, random_state=i)
    clf = SVC(kernel="linear", C=1)
    clf.fit(kma_train, kma_y_train)
    scores.append(clf.score(kma_test,kma_y_test))


scores = np.array(scores)
values, counts = np.unique(scores, return_counts=True)

fig = plt.figure(frameon=False)
ax = fig.add_subplot(111)
ax.bar(values, counts, 0.05)
ax.set_title(f"Mean: {np.mean(scores):.2f} σ: {np.std(scores):.2f}")

ax.set_xlabel("score")
ax.set_ylabel("count")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels)
fig.savefig("plots/kma_only_scores.pdf")











