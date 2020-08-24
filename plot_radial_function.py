import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

grid = np.loadtxt("dat/grid.dat")

plt.figure(1)
# for i in range(7):
for i in range(1):
    fname = "dat/eigenvector_function_%i.dat" % i
    eigenfunction = np.loadtxt(fname)

    x = grid
    f = eigenfunction[:, 0]
    f[1:] = f[1:]/x[1:]
    f2 = f*f
    normed_f = f2/sum(f2)
    plt.plot(x, normed_f)



plt.xlim([0.0, 20.0])
plt.show()