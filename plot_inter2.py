import numpy as np
import pandas as pd
from sklearn.datasets import load_iris
from ipywidgets import interact, Select
import matplotlib.pyplot as plt


def show_plot(col1, col2):
    plt.figure(figsize=(4.5, 4.5))
    plt.scatter(data[col1], data[col2], c=iris.target)
    plt.xlabel('')
    plt.ylabel('')
    plt.show()


iris = load_iris()
data = pd.DataFrame(iris.data, columns=iris.feature_names)

w1 = Select(description='X-axis:', options=data.columns, rows=4,)
w2 = Select(description='Y-axis:', options=data.columns, rows=4,)
interact(show_plot, col1=w1, col2=w2)
