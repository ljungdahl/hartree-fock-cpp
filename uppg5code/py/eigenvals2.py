from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from Hmatrix import Hmat,H_GQint
from Bmatrix import Bmat,B_GQint
from globals import * #imports l,Z,N,k,rmin,rmax - the global parameters.

Hh = Hmat(N,k,rmin,rmax,l,Z)
Bb = Bmat(N,k,rmin,rmax)

#print Hh[-1][-1], Bb[-1][-1] #SEEMS TO WORK

w, vr = LA.eig(Hh,Bb)
tt = knotSeq(k,N,rmin,rmax)
xx = linspace(rmin,rmax,1000)
minidx = argmin(w)

print w

"""
for i in arange(1,len(w)+1):
    plt.plot(xx,B(k,tt,i+1,xx,0))

plt.show()
"""


basisMatrixn1 = zeros((len(xx),len(w)))
minCoeffs = vr[:,minidx]
#print minCoeffs
for i in arange(0,len(xx)):
    for j in arange(0,len(w)):
        basisMatrixn1[i][j] = minCoeffs[j]*B(k,tt,j+2,xx[i],0)

superposn1 = sum(basisMatrixn1,axis=1)
superpos2n1 = pow(abs(superposn1),2)
normsp2n1 = superpos2n1/sum(superpos2n1)

w2 = delete(w,minidx)
minidx2 = argmin(w2)

minindices = w.argsort()[:3]
print minindices
print minidx

basisMatrixn2 = zeros((len(xx),len(w)))
minCoeffs = vr[:,minidx2]
#print minCoeffs
for i in arange(0,len(xx)):
    for j in arange(0,len(w)):
        basisMatrixn2[i][j] = minCoeffs[j]*B(k,tt,j+2,xx[i],0)

superposn2 = sum(basisMatrixn2,axis=1)
superpos2n2 = pow(abs(superposn2),2)
normsp2n2 = superpos2n1/sum(superpos2n2)

a0 =1

plt.plot(xx,normsp2n1,label='n=1, l=0, E='+str(w[minidx]),linewidth=2)
plt.plot(xx,normsp2n1,label='n=1, l=0, E='+str(w2[minidx2]),linewidth=2)
plt.axvline(a0,0,10,color='r',linestyle='--',linewidth=2)
plt.legend()
plt.show()

