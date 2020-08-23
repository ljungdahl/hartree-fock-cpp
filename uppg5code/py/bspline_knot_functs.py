from numpy import *
from scipy.interpolate import BSpline

#BSpline function
def B(order,knots,number,space,deriv):
    c = zeros(len(knots))
    c[number-1] = 1
    Bspl = BSpline(knots,c,order-1,extrapolate=False)
    if deriv == 0:
        rBspl = Bspl(space)
    else:
        Bspl = Bspl.derivative(deriv)
        rBspl = Bspl(space)
    
    nanIdx = isnan(rBspl)
    rBspl[nanIdx] = 0
    return rBspl

#knotpoint construction (linearly spaced!)
def knotSeq(order,totalknots,RMIN,RMAX):
    btweenpts = totalknots-(2*(order-1))
    RMINARR = RMIN*ones(order-1)
    RMAXARR = RMAX*ones(order-1)
    RBTWN = linspace(RMIN,RMAX,btweenpts)
    #RBTWN1 = linspace(RMIN,RMAX/5.0,btweenpts/2)
    #RBTWN2 = linspace(RMAX-(4.0/5.0)*RMAX,RMAX,(btweenpts/2)+1) #The sequence between RMIN and RMAX
    t = append(RMINARR,RBTWN)
    t = append(t, RMAXARR)
    #t = append(RMINARR, RBTWN1)
    #t = append(t,RBTWN2[1:])
    #t = append(t, RMAXARR)
    
    return t

