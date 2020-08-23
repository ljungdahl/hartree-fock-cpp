from numpy import *
import matplotlib.pyplot as plt

def aPsi(x,n,l):
    n1l0 = exp(-1.0*x)#(1.0/sqrt(pi))*exp(-1.0*x)
    n2l0 = (2-x)*exp(-0.5*x)#(1.0/(4*sqrt(2*pi)))*(2-x)*exp(-0.5*x)
    n3l0 = (27-18*x+2*pow(x,2))*exp(-(1.0/3.0)*x)#(1.0/(81*sqrt(3*pi)))*(27-18*x+2*pow(x,2))*exp(-(1.0/3.0)*x)
    
    n2l1 = x*exp(-0.5*x)
    n3l1 = (6-x)*x*exp(-(1.0/3.0)*x)

    n3l2 = pow(x,2)*exp(-(1.0/3.0)*x)

    if l == 0:
        if n == 1:
            return n1l0
        elif n == 2:
            return n2l0
        elif n == 3:
            return n3l0
    elif l == 1:
        if n == 1:
            return ones(len(x))
        elif n == 2:
            return n2l1
        elif n == 3:
            return n3l1
    elif l == 2 and n == 3:
        return n3l2
    else:
        return ones(len(x))
    
"""
xx = linspace(0,10,1000)

plt.plot(xx,pow(xx,2)*pow(aPsi(xx,1),2))
plt.plot(xx,pow(xx,2)*pow(aPsi(xx,2),2))
plt.show()
"""
