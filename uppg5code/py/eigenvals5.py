from numpy import *
from cycler import cycler
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from Hmatrix import Hmat,H_GQint
from Bmatrix import Bmat,B_GQint
from globals import * #imports l,Z,N,k,rmin,rmax - the global parameters.
from analyticpsi import aPsi



plt.rc('text', usetex=True)
plt.rc('font', family='serif')
params = {'text.latex.unicode' : True, 'font.size' : 20}
plt.rcParams.update(params)

plt.rc('lines', linewidth=4)
plt.rc('axes', prop_cycle=(cycler('color', ['c', 'r', 'm','y','b']) +
                           cycler('linestyle', ['--','--','--','--','--'])))

HtGS = (13.60000000/27.211400000000) # ~= 0.5 Hartree

Hh = Hmat(N,k,rmin,rmax,l,Z)
Bb = Bmat(N,k,rmin,rmax)

#print Hh[-1][-1], Bb[-1][-1] #SEEMS TO WORK

w, vr = LA.eig(Hh,Bb)
tt = knotSeq(k,N,rmin,rmax)
xx = linspace(rmin,rmax,1000)
minidx = argwhere(w < 0)
minidx = minidx[:,0]
#print minidx, len(minidx)

#print w
aPsi2n1 = pow(abs(xx*aPsi(xx,1,0)),2)
aPsi2n1 = aPsi2n1/sum(aPsi2n1)
aPsi2n2 = pow(abs(xx*aPsi(xx,2,0)),2)
aPsi2n2 = aPsi2n2/sum(aPsi2n2)
aPsi2n3 = pow(abs(xx*aPsi(xx,3,0)),2)
aPsi2n3 = aPsi2n3/sum(aPsi2n3)

if Z == 1:
    plt.plot(xx,aPsi2n1,linestyle='-',color='k',linewidth=1)
    plt.plot(xx,aPsi2n2,linestyle='-',color='k',linewidth=1)
    if len(minidx) > 2:
        plt.plot(xx,aPsi2n3,linestyle='-',color='k',linewidth=1)

"""
for i in arange(1,len(w)+1):
    plt.plot(xx,B(k,tt,i+1,xx,0))

plt.show()
"""
basisMatrix = zeros((len(xx),len(w)))
for i in arange(0,len(xx)):
    for j in arange(0,len(w)):
        basisMatrix[i][j] = B(k,tt,j+2,xx[i],0)

nn = 1
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
    calclabelstr = r'$E_{%i_c} = %.4f$, ' % (nn,eigenvalue)
    anaeig = (-HtGS*Z**2)/nn**2
    analabelstr = r'$E_{%i_a} = %.4f$' % (nn,anaeig)
    labelstr = calclabelstr+analabelstr
    plt.plot(xx,normsp2,label=labelstr)
    nn = nn+1
#plt.axvline(1,0,10,color='r',linestyle='--',linewidth=2)
plt.title(r'The bound states for $N=%i,k=%i,Z=%i,l=%i,R_{max}=%.1f$' % (N,k,Z,l,rmax))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$r/a_0$')
plt.ylabel(r'$r^2 |\Psi|^2$')
plt.grid()
plt.ylim
plt.legend()
plt.show()

