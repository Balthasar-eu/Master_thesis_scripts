#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 05:29:28 2019

@author: balthasar
"""

import numpy as np

import autosklearn.classification
import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics

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
"ERR1100976" : 16
}


csv = np.loadtxt("all_genomes.csv", delimiter=",")

print(len(csv[0]))

csv = csv[:,np.std(csv, axis=0) != 0]

#print(csv)

print(len(csv[0]))

newv = {v: (0 if virulence[v]<50 else 1) for v in virulence}

with open("csvs/keyfile.txt","r") as f:
    keys = [line for line in f]

Y = np.array([newv[k[:-2]] for k in keys], dtype=int)

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(csv, Y, random_state=1)
automl = autosklearn.classification.AutoSklearnClassifier()
automl.fit(X_train, y_train)
y_hat = automl.predict(X_test)
print("Accuracy score", sklearn.metrics.accuracy_score(y_test, y_hat))






