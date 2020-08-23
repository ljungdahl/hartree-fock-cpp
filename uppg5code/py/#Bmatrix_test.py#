from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from globals import * #imports l,Z,N,k,rmin,rmax - the global parameters.

def B_GQint(k,t,i,j):
    if k == 4:
    #n = 4 weights and abscissae (Source https://pomax.github.io/bezierinfo/legendre-gauss.html)
        weights = array([0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538]) 
        absc = array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
    elif k == 6:
    #n = 6 weights and abscissae (4.4*1e-6 error for single Bspline)
        weights = array([0.3607615730481386,0.3607615730481386,0.4679139345726910,0.4679139345726910,0.1713244923791704,0.1713244923791704])
        absc = array([0.6612093864662645,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,-0.932469514203152,0.932469514203152])
    
    I1 = 0
    
    #make an interval for integration between relevant knotpts
    for idx1 in arange(max(i,j),min(i,j)+k):
        
        
        a = t[idx1-1]
        b = t[idx1]
        #construct evaluation points on the interval
        x_gq = zeros(len(absc))
        hej = 1
        if abs(b-a) > 0:
        #if hej == 1:
        #perform integration on interval
            for idx2 in arange(0,len(x_gq)):
                zz = absc[idx2]
                x_gq[idx2] = (zz*(b-a)*0.5)+(a+b)*0.5 #We transform the zz
                                                      #in [-1,1] to x_gq in [a,b].

            ff = B(k,t,i,x_gq,0)*B(k,t,j,x_gq,0)
            gq_int = (b-a)*0.5*dot(weights,ff) #the integral for a given interval
        else:
            gq_int = 0
            
        I1 = I1+gq_int #summing up the integrals over the intervals
    
    
    return I1
#End of H_GQint()




#knotpoints 




#Construct B.
def Bmat(N,k,rmin,rmax):
    dim = N-k-2
    t = knotSeq(k,N,rmin,rmax)
    BB = zeros((dim,dim))

    for ii in range(0,dim):  
        for jj in range(0,dim):
            i = ii+2
            j = jj+2
            if abs(i-j) <= k-1:
                BB[ii][jj] = B_GQint(k,t,i,j)
            else:
                BB[ii][jj] = 0.0

    return BB


