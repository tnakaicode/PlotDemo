#!/usr/bin/env python
# coding: utf-8
# https://gist.github.com/dcf98f99479e18d7a9e0642600d2893c.git
# # Scatterplot with histogram of differences


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator

from base import plot2d

obj = plot2d()

n = 100
rng = np.random.default_rng()
x = rng.uniform(0.3, 0.85, n)
y = x - rng.normal(.05, .05, n)

lim = .3, .9
obj.axs.set(xlim=lim, ylim=lim)
# Forthcoming ax.axline makes this easier
obj.axs.plot(lim, lim, c=".7", dashes=(4, 2), zorder=0)
obj.axs.scatter(x, y)
obj.axs.set(xlabel="x", ylabel="y")

obj.fig.subplots_adjust(right=.7, top=.7)
obj.axs.spines["right"].set_visible(False)
obj.axs.spines["top"].set_visible(False)

hist, bins = np.histogram(x - y)
w = bins[1] - bins[0]

plot_extents = -.2, .2, 0, hist.max()

transform = Affine2D().scale(60, 1).rotate_deg(-45)
helper = floating_axes.GridHelperCurveLinear(
    transform, plot_extents, grid_locator1=MaxNLocator(4))
inset = floating_axes.FloatingAxes(
    obj.fig, [.65, .65, .35, .35], grid_helper=helper)

bar_ax = inset.get_aux_axes(transform)
bar_ax.bar(bins[:-1], hist, w)
bar_ax.plot((0, 0), (0, hist.max()), c=".7", dashes=(4, 2))

inset.axis["left"].set_visible(False)
inset.axis["right"].set_visible(False)
inset.axis["top"].set_visible(False)

axis = inset.axis["bottom"]
axis.major_ticks.set_tick_out(True)

obj.fig.add_axes(inset)
#obj.fig.set_size_inches(18.5, 10.5)
obj.fig.set_size_inches(6.0, 6.0)
obj.SavePng()
#f.savefig("difference_hist.png", dpi=200)
