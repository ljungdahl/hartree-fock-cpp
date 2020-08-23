from numpy import *
from cycler import cycler
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from Hmatrix import Hmat,H_GQint
from Bmatrix import Bmat,B_GQint
from globals import * #imports l,Z,N,k,rmin,rmax,lcompflag - the global parameters.
from analyticpsi import aPsi



plt.rc('text', usetex=True)
plt.rc('font', family='serif')
params = {'legend.fontsize': 15, 'text.latex.unicode' : True, 'font.size' : 20}
plt.rcParams.update(params)

plt.rc('lines', linewidth=2)
plt.rc('axes', prop_cycle=(cycler('color', ['c', 'r', 'm','y','b']) +
                           cycler('linestyle', ['-','-','-','-','-'])))

HtGS = (13.6/27.211385) # ~= 0.5 Hartree

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
aPsi2n1 = pow(abs(xx*aPsi(xx,1,l)),2)
aPsi2n1 = aPsi2n1/sum(aPsi2n1)
aPsi2n2 = pow(abs(xx*aPsi(xx,2,l)),2)
aPsi2n2 = aPsi2n2/sum(aPsi2n2)
aPsi2n3 = pow(abs(xx*aPsi(xx,3,l)),2)
aPsi2n3 = aPsi2n3/sum(aPsi2n3)

apsilab = r'Analytic $r^2 |R_{n,l}(r)|^2$'
if Z == 1:
    if l == 0:
        plt.plot(xx,aPsi2n1,linestyle='--',color='k',linewidth=4,label=apsilab)
        plt.plot(xx,aPsi2n2,linestyle='--',color='k',linewidth=4) #label=apsilab)
        if len(minidx) > 2:
            plt.plot(xx,aPsi2n3,linestyle='--',color='k',linewidth=4) #label=apsilab)
    elif l == 1:
        plt.plot(xx,aPsi2n2,linestyle='--',color='k',linewidth=4) #label=apsilab)
        plt.plot(xx,aPsi2n3,linestyle='--',color='k',linewidth=4) #label=apsilab)
    elif l == 2:
        plt.plot(xx,aPsi2n3,linestyle='--',color='k',linewidth=4) #label=apsilab)

"""
for i in arange(1,len(w)+1):
    plt.plot(xx,B(k,tt,i+1,xx,0))

plt.show()
"""
basisMatrix = zeros((len(xx),len(w)))
for i in arange(0,len(xx)):
    for j in arange(0,len(w)):
        basisMatrix[i][j] = B(k,tt,j+2,xx[i],0)

if l == 0:
    nn = 1
elif l == 1:
    nn = 2
elif l == 2:
    nn = 3
    
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
    if nn < 4:
        plt.plot(xx,normsp2,label=labelstr)
    nn = nn+1

fname = 'plot_N%i_k%i_l%i_Z%i_rmax%2.0f.eps' % (N,k,l,Z,rmax)
if lcompflag == 1:
    fname = 'lcomp_'+fname

path = '/home/anton/s/comp_phys/a5/tex/'
fullfname = path+fname
#plt.axvline(1,0,10,color='r',linestyle='--',linewidth=2)
plt.title(r'The bound states for $Z=%i,l=%i,R_{max}=%.1f$' % (Z,l,rmax))
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$r/a_0$')
plt.ylabel(r'$r^2 |R_{n l}|^2$')
if rmax < 16:
    plt.xticks(arange(rmin, rmax+1.0,1.0))
plt.grid()
plt.ylim
plt.legend()
plt.savefig(fullfname,dpi=1000, bbox_inches='tight')

print 'Done! Saved to %s ' % (fullfname)
plt.show()

