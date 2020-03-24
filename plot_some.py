import matplotlib.pyplot as plt
import numpy as np

from base import plot2d

num = 10 * 10**6

px = np.linspace(0, 1, num) * 1000

obj = plot2d(aspect="auto")
obj.axs.plot(px, 0.1 * px, label="0.1")
obj.axs.plot(px, 0.2 * px, label="0.2")
obj.axs.plot(px, 0.3 * px, label="0.3")

# 
# UserWarning: Creating legend with loc="best" can be slow with large amounts of data. 
# self.fig.savefig(pngname)
# 

#obj.axs.plot(px, 0.4 * px, label="0.4")
#obj.axs.plot(px, 0.5 * px, label="0.5")
#obj.axs.plot(px, 0.6 * px, label="0.6")
obj.axs.legend()
obj.SavePng_Serial()
print("ok")
