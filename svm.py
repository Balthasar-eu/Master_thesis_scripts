#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:44:33 2019

@author: balthasar
"""

from sklearn.svm import SVC
import numpy as np
import matplotlib.pyplot as plt
import sklearn
from sklearn.model_selection import GridSearchCV
from sklearn.decomposition import PCA


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

parameters = [
        {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}
        ]

csv = np.loadtxt("data/all_genomes.csv", delimiter=",")
arla_csv = np.loadtxt("data/all_arla_genomes.csv", delimiter=",")

len_before_trimming = len(csv[0])
std_genes = np.std(csv, axis=0)
print(std_genes)

print(csv.shape)
print(arla_csv.shape)

std_sorted = np.sort(std_genes)
std_sumcount = np.arange(1,len(std_sorted)+1)

arla_csv = arla_csv[:,np.std(csv, axis=0) > 15]
csv = csv[:,np.std(csv, axis=0) > 15]

pca = PCA(n_components=3)
pca.fit(csv.T)

print(pca.explained_variance_ratio_)
print(pca.components_)
print(pca.components_.shape)

np.savetxt("data/pca_components.csv", pca.components_, delimiter=",")

pca = PCA(n_components=3)
pca.fit(np.vstack((csv,arla_csv)).T)

np.savetxt("data/pca_components_arla.csv", pca.components_, delimiter=",")

len_after_trimming = len(csv[0])

newv = {sn: (1 if v<25 else (2 if v < 50 else (3 if v < 60 else 4)) ) for sn,v in virulence.items()}

count = {1:0,2:0,3:0,4:0}
for k in newv:
    count[newv[k]] += 1

with open("data/keyfile.txt","r") as f:
    keys = [line for line in f]

Y = np.array([newv[k[:-10]] for k in keys], dtype=int)

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=1)

clf = GridSearchCV(SVC(), parameters, cv=5)

#clf = SVC(gamma='auto')
clf.fit(X_train, y_train)
print(clf.score(X_test,y_test))
with open("data/svm_results.txt","w") as f:
    f.write("Before trimming: "+ str(len_before_trimming)+"\n")
    f.write("After trimming: "+ str(len_after_trimming)+"\n")
    f.write("Count of classes:\n" + str(count)+"\n")
    f.write("Train length: "+ str(len(X_train))+"\n")
    f.write("Test length: " + str(len(X_test))+"\n")
    f.write("Accuracy: "+ str(clf.score(X_test,y_test)) +"\n")
    f.write("Parameters:" + str(clf.best_params_)+"\n")
    f.write("Gridsearch results: " + str(clf.cv_results_)+"\n")
    f.write(str(clf.predict(X_test))+"\n")
    f.write(str(y_test) +"\n" )
    f.write(str(clf.predict(arla_csv)) + "\n")
#    f.write(str(clf.coef_))
#    f.write(str(clf.intercept_))

scores = []

for i in range(100):
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=i)
    clf = SVC(kernel="linear", C=1)
    clf.fit(X_train, y_train)
    scores.append(clf.score(X_test,y_test))
    
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
fig.savefig("plots/scores.pdf")

