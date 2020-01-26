#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 12:08:01 2019

@author: balthasar
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA

plt.close("all")

###############################################################################
################################## 3D #########################################

# we create 40 separable points
np.random.seed(0)
X = np.r_[np.random.randn(20, 3) - [2, 2, 2], np.random.randn(20, 3) + [2, 2, 2]]
Y = [0] * 20 + [1] * 20

# fit the model
clf = svm.SVC(kernel='linear')
clf.fit(X, Y)

# get the separating hyperplane

#n = np.cross(np.array([1,1,1]),clf.coef_[0])

n = clf.coef_[0]
d = clf.intercept_
xx, yy = np.meshgrid([-10,10], [-10,10])
zz = (-n[0] * xx - n[1] * yy - d) * 1/n[2]
# plot the line, the points, and the nearest vectors to the plane
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=Y, cmap=plt.cm.Paired)
ax.axis('tight')

fig.savefig("plots/testdata.svg")

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=Y, cmap=plt.cm.Paired)
ax.axis('tight')

ax.plot_surface(xx,yy, zz, alpha=0.2)
fig.savefig("plots/testdata_hyperplane.svg")

pca = PCA()

Xnew = pca.fit_transform(X)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

transformed_stacked_mesh = pca.transform(np.dstack([np.concatenate([row for row in dim]) for dim in [xx,yy,zz]])[0])
split_tsm = np.vsplit(transformed_stacked_mesh, 2)
meshlist = [np.vstack([s[:,i] for s in split_tsm]) for i in range(3)]
ax.plot_surface(meshlist[0],meshlist[1],meshlist[2], alpha=0.2)
ax.scatter(Xnew[:, 0], Xnew[:, 1], Xnew[:, 2], c=Y, cmap=plt.cm.Paired)
ax.axis('tight')

fig.savefig("plots/testdata_pca_hyperplane.svg")



plt.close("all")