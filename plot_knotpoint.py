import numpy as np
import matplotlib.pyplot as plt


grid = np.loadtxt("grid.dat")
knotpts = np.loadtxt("knotpoints.dat")

plt.figure(1)
plt.plot(grid, np.zeros(len(grid)))
plt.plot(knotpts, np.zeros(len(knotpts)), 'kd')
plt.show()


