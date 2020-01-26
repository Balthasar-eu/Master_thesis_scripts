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
"10CEB535LM" : 7,
"10CEB540LM" : 7,
"10CEB550LM" : 50.9,
"10CEB552LM" : 50.9,
"10CEB553LM" : 45,
"10CEB554LM" : 45,
"10CEB559LM" : 33.3,
"10CEB560LM" : 33.3,

}

parameters = [
        {'kernel': ['rbf'], 'gamma': np.logspace(-5,5,11), 'C': [1, 10, 100, 1000]},
        {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}
        ]

csv = np.loadtxt("supplemental_files/all_genomes.csv", delimiter=",")

len_before_trimming = len(csv[0])
csv[csv<50] = 0
std_genes = np.std(csv, axis=0)
print(std_genes)

std_sorted = np.sort(std_genes)
std_sumcount = np.arange(1,len(std_sorted)+1)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(std_sorted, std_sumcount, label='cumulative std count')
fig.savefig("results/std_count.pdf")

csv = csv[:,np.std(csv, axis=0) > 15]

len_after_trimming = len(csv[0])

newv = {v: (0 if virulence[v]<25 else (1 if virulence[v] < 50 else (2 if virulence[v] < 60 else 3)) ) for v in virulence}

count = {0:0,1:0,2:0,3:0}
for k in newv:
    count[newv[k]] += 1

with open("supplemental_files/keyfile.txt","r") as f:
    keys = [line for line in f]

Y = np.array([newv[k[:-10]] for k in keys], dtype=int)

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=1)

clf = GridSearchCV(SVC(), parameters, cv=5)

#clf = SVC(gamma='auto')
clf.fit(X_train, y_train)
print(clf.score(X_test,y_test))
with open("results/svm_manually/results_4_classes.txt","w") as f:
    f.write("Before trimming: "+str(len_before_trimming)+"\n")
    f.write("After trimming: "+str(len_after_trimming)+"\n")
    f.write("Count of classes:\n" + str(count)+"\n")
    f.write("Train length: "+ str(len(X_train))+"\n")
    f.write("Test length: " + str(len(X_test))+"\n")
    f.write("Accuracy: "+ str(clf.score(X_test,y_test)) +"\n")
    f.write("Parameters:" + str(clf.best_params_)+"\n")
    f.write("Gridsearch results: " + str(clf.cv_results_)+"\n")

################################################ 2 classes

newv = {v: (0 if virulence[v]<50 else 1) for v in virulence}

count = {0:0,1:0,2:0,3:0}
for k in newv:
    count[newv[k]] += 1

with open("supplemental_files/keyfile.txt","r") as f:
    keys = [line for line in f]

Y = np.array([newv[k[:-10]] for k in keys], dtype=int)

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=1)

clf = GridSearchCV(SVC(), parameters, cv=5)

#clf = SVC(gamma='auto')
clf.fit(X_train, y_train)
print(clf.score(X_test,y_test))
with open("results/svm_manually/results_2_classes.txt","w") as f:
    f.write("Before trimming: "+str(len_before_trimming)+"\n")
    f.write("After trimming: "+str(len_after_trimming)+"\n")
    f.write("Count of classes:\n" + str(count)+"\n")
    f.write("Train length: "+ str(len(X_train)) + "\n")
    f.write("Test length: " + str(len(X_test)) + "\n")
    f.write("Accuracy: "+ str(clf.score(X_test,y_test)) +"\n")
    f.write("Parameters:" + str(clf.best_params_) + "\n")
    f.write("Gridsearch results: " + str(clf.cv_results_) +"\n")











