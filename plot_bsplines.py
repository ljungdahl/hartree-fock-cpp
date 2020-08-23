import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

def B(x, k, i, t):
    if k == 0:
        return 1.0 if t[i] <= x < t[i+1] else 0.0
    if t[i+k] == t[i]:
        c1 = 0.0
    else:
        c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
    if t[i+k+1] == t[i+1]:
        c2 = 0.0
    else:
        c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
    return c1 + c2

def bspline(x, t, c, k):
    n = len(t) - k - 1
    assert (n >= k+1) and (len(c) >= n)
    return sum(c[i] * B(x, k, i, t) for i in range(n))


data = np.loadtxt("bsplines.dat")
data2 = np.loadtxt("bsplines_derivatives.dat")
knotpts = np.loadtxt("knotpoints.dat")

plotDerivs = True


x = data[:, 0]
x2 = data2[:, 0]
bsplines = data[:, 1:]
dBsplines = data2[:, 1:]


order = 4
k = order

numKnotPts = len(knotpts) # n+k+1
numPhysPts = numKnotPts - 2 * (k-1)
numBsplines = numKnotPts-k #n+k+1-k = n+1
#t = knotpts[k-1:len(knotpts)-(k-1)]
t = knotpts

print(t)
for bsplIdx in range(numBsplines):
    c = np.zeros(numBsplines)
    c[bsplIdx] = 1.0
    B = spi.BSpline(t, c, k, False)
    dB = B.derivative()
    plt.plot(x, B(x))
    plt.plot(x, dB(x), 'b--')

# def getBasisDerivative(knotPts, i, order, xx):
#     k = order
#     t = knotPts
#     B_i_k_1 = spi.BSpline.basis_element(knotPts[i:i+k+1-1], True)
#     B_i_p1_km1 = spi.BSpline.basis_element(knotPts[i+1:i+1+k+1-1], True)
#     dB_i_k = (k-1)*(B_i_k_1(xx)/(t[i+k-1]-t[i])-B_i_p1_km1(xx)/(t[i+k]-t[i+1]))
#     return dB_i_k

# ## Basis elements nice! Can I get the same with regular so I can get derivative?
# pythonBsplines = []
# plt.figure(1)
# for bsplIdx in range(numBsplines):
#     t = knotpts[bsplIdx:bsplIdx+k+1]
#     B = spi.BSpline.basis_element(t, False)
#     pythonBsplines.append(B)
#     dB = getBasisDerivative(knotpts, bsplIdx, k, x)
#     plt.plot(x, dB, 'b--')
#     plt.plot(x, B(x))

plt.figure(2)
# myBsplIdx = 1
# plt.plot(x, bsplines[:, myBsplIdx])
# plt.plot(x2, dBsplines[:, myBsplIdx], 'b--')
# myBsplIdx = numBsplines-1-1
# plt.plot(x, bsplines[:, myBsplIdx])
# plt.plot(x2, dBsplines[:, myBsplIdx], 'b--')
for i in range(bsplines.shape[1]):
    #print(i)
    plt.plot(x, bsplines[:, i])

    if(plotDerivs):
        plt.plot(x2, dBsplines[:, i], 'b--')


knotptsConst = 0.01*np.ones(len(knotpts))
plt.plot(knotpts, knotptsConst, 'kd')
# plt.plot(x, B(x), 'b--')
# plt.plot(x, pythonBsplines(x), 'b--')
# plt.plot(x, bsplinevals, 'b--')

plt.show()


