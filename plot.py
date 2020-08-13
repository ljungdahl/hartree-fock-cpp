import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("bsplines.dat")

plt.figure(1)
x = data[:, 0]
bsplines = data[:, 1:]

# plt.plot(x, y)

for i in range(bsplines.shape[1]):
    print(i)
    plt.plot(x, bsplines[:, i])

# plt.plot(data)
plt.show()
