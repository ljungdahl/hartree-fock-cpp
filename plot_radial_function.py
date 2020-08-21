import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

data = np.loadtxt("bsplines.dat")

radfunc_data = np.loadtxt("radialfunction.dat")

x = data[:, 0]

f = radfunc_data/x
f2 = f*f

plt.plot(x, f2)

plt.show()


