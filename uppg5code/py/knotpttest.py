from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from globals import * #imports l,Z,N,k,rmin,rmax - the global parameters.


#n = 6 weights and abscissae (4.4*1e-6 error for single Bspline)
#weights = array([0.3607615730481386,0.3607615730481386,0.4679139345726910,0.4679139345726910,0.1713244923791704,0.1713244923791704])
#absc = array([0.6612093864662645,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,-0.932469514203152,0.932469514203152])
t = knotSeq(k,N,rmin,rmax)
print t[0:10]
xx = linspace(rmin,rmax,1000)
#n = 4 weights and abscissae (Source https://pomax.github.io/bezierinfo/legendre-gauss.html)
weights = array([0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538]) 
absc = array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
i = 2
j = 3

I1 = 0
for idx1 in arange(max(i,j)+(k-2),min(i,j)+(k-2)+k-1):
    a = t[idx1]
    b = t[idx1+1]
    print a,b
    #construct evaluation points on the interval
    x_gq = zeros(len(absc))
    #perform integration on interval
    for idx2 in arange(0,len(x_gq)):
        zz = absc[idx2]
        x_gq[idx2] = (zz*(b-a)*0.5)+(a+b)*0.5 #We transform the zz in [-1,1] to x_gq in [a,b].

    print x_gq
    ff = B(k,t,i,x_gq,0)*B(k,t,j,x_gq,0)
    gq_int = (b-a)*0.5*dot(weights,ff) #the integral for a given interval
        
        
    I1 = I1+gq_int #summing up the integrals over the intervals

print I1
plt.plot(xx,B(k,t,1,xx,0))
plt.plot(xx,B(k,t,i,xx,0))
plt.plot(xx,B(k,t,j,xx,0))
plt.xlim([0,0.5])
plt.show()




## Test that compares with simpson integration, it works!!!!
"""
tmin = t[max(i,j)]
tmax = t[min(i,j)+k-1]
print tmin,tmax
nn = 10000
h = abs(tmax-tmin)/nn
#print h
simp_integspace = linspace(tmin,tmax,nn)
xx_s = simp_integspace
I2 = 0
for idx3 in arange(0,len(simp_integspace)):
    ffx = B(k,t,i,xx_s[idx3],0)*B(k,t,j,xx_s[idx3],0)
    I2 = I2+ffx*h

print I2
"""
#plt.plot(xx,B(k,t,i,xx,0))
#plt.plot(xx,B(k,t,j,xx,0))
#plt.show()

