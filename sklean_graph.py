import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from sklearn.datasets import load_iris
from sklearn import tree
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = plot2d(aspect="equal")

    clf = tree.DecisionTreeClassifier(random_state=0)
    iris = load_iris()
    clf = clf.fit(iris.data, iris.target)
    tree.plot_tree(clf)
    obj.SavePng()
