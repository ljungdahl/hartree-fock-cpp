import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

grid = np.loadtxt("dat/grid.dat")
knots = np.loadtxt("dat/knotpoints.dat")

# B5 = np.loadtxt("dat/B_5.dat")
# dB5 = np.loadtxt("dat/dB_5.dat")
# d2B5 = np.loadtxt("dat/dB2_5.dat")

B0 = np.loadtxt("dat/B_0.dat")
dB0 = np.loadtxt("dat/dB_0.dat")
d2B0 = np.loadtxt("dat/dB2_0.dat")

B2 = np.loadtxt("dat/B_2.dat")
dB2 = np.loadtxt("dat/dB_2.dat")
d2B2 = np.loadtxt("dat/dB2_2.dat")

B5 = np.loadtxt("dat/B_5.dat")
dB5 = np.loadtxt("dat/dB_5.dat")
d2B5= np.loadtxt("dat/dB2_5.dat")

B7 = np.loadtxt("dat/B_7.dat")
dB7 = np.loadtxt("dat/dB_7.dat")
d2B7= np.loadtxt("dat/dB2_7.dat")

B11 = np.loadtxt("dat/B_11.dat")
dB11 = np.loadtxt("dat/dB_11.dat")
d2B11= np.loadtxt("dat/dB2_11.dat")

Bsplines = []
Bsplines.append(B11)
Bsplines.append(dB11)
Bsplines.append(d2B11)

plt.figure(1)

for bspline in Bsplines:
    plt.plot(grid, bspline)
    plt.plot(knots, np.zeros(len(knots)), 'kd')

plt.show()