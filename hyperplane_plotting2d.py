#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 12:08:01 2019

@author: balthasar
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm

plt.close("all")

###############################################################################
################################## 2D #########################################
# we create 40 separable points
np.random.seed(10)
X = np.r_[np.random.randn(20, 2) - [2, 2], np.random.randn(20, 2) + [2, 2]]
Y = [0] * 20 + [1] * 20

# fit the model
clf = svm.SVC(kernel='linear')
clf.fit(X, Y)

# get the separating hyperplane
w = clf.coef_[0]
a = -w[0] / w[1]
xx = np.linspace(-5, 5)
yy = a * xx - (clf.intercept_[0]) / w[1]

# plot the parallels to the separating hyperplane that pass through the
# support vectors
b = clf.support_vectors_[0]
yy_down = a * xx + (b[1] - a * b[0])
b = clf.support_vectors_[-1]
yy_up = a * xx + (b[1] - a * b[0])

# plot the line, the points, and the nearest vectors to the plane
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xticks([])
ax.set_yticks([])

ax.set_ylim([-5,5])
ax.set_xlim([-5,5])

ax.scatter(X[:, 0], X[:, 1], c=Y, cmap=plt.cm.Paired)
fig.savefig("plots/2dtestdata.pdf")


ax.plot(xx, yy, 'k-')
ax.plot(xx, yy_down, 'k--')
ax.plot(xx, yy_up, 'k--')
ax.scatter(clf.support_vectors_[:, 0], clf.support_vectors_[:, 1],
            s=80, facecolors='none')
fig.savefig("plots/2dtestdata_hyperplane.pdf")
#plt.close("all")

