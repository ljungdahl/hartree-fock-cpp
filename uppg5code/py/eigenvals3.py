from numpy import *
from cycler import cycler
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from Hmatrix import Hmat,H_GQint
from Bmatrix import Bmat,B_GQint
from globals import * #imports l,Z,N,k,rmin,rmax - the global parameters.

plt.rc('lines', linewidth=3)
plt.rc('axes', prop_cycle=(cycler('color', ['b', 'r', 'k']) +
                           cycler('linestyle', ['--','--','--'])))

Hh = Hmat(N,k,rmin,rmax,l,Z)
Bb = Bmat(N,k,rmin,rmax)

#print Hh[-1][-1], Bb[-1][-1] #SEEMS TO WORK

w, vr = LA.eig(Hh,Bb)
tt = knotSeq(k,N,rmin,rmax)
xx = linspace(rmin,rmax,1000)
minidx = w.argsort()[:3]

print w

"""
for i in arange(1,len(w)+1):
    plt.plot(xx,B(k,tt,i+1,xx,0))

plt.show()
"""
basisMatrix = zeros((len(xx),len(w)))
for i in arange(0,len(xx)):
    for j in arange(0,len(w)):
        basisMatrix[i][j] = B(k,tt,j+2,xx[i],0)

for mm in arange(0,len(minidx)):

    #for i in arange(0,len(xx)):
    #    for j in arange(0,len(w)):
    #        basisMatrix[i][j] = (k,tt,j+2,xx[i],0)

    coeffdot = vr[:,minidx[mm]]*basisMatrix
    #print coeffdot
    superpos = sum(coeffdot,axis=1)
    #print superpos
    spsqrd = pow(abs(superpos),2)
    normsp2 = spsqrd/sum(spsqrd)
    eigenvalue = float(real(w[minidx[mm]]))
    print eigenvalue
    labelstr = 'E = {0:.2f}'.format(eigenvalue)
    plt.plot(xx,normsp2,label=labelstr)

    
#plt.axvline(1,0,10,color='r',linestyle='--',linewidth=2)
plt.grid()
plt.legend()
plt.show()

