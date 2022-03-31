import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from sklearn.datasets import load_iris
from sklearn import tree
from linecache import getline, clearcache
import argparse

sys.path.append(os.path.join("./"))
from base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = plot2d(aspect="equal")

    clf = tree.DecisionTreeClassifier(random_state=0)
    iris = load_iris()
    clf = clf.fit(iris.data, iris.target)
    tree.plot_tree(clf)
    obj.SavePng()
