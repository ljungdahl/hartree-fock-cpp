from numpy import *
import scipy.interpolate as si
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
import scipy.linalg as LA

#BSpline function
def B(order,knots,number,space):
    c = zeros(len(knots))
    c[number-1] = 1
    Bspl = BSpline(knots,c,order,extrapolate=False)
    rBspl = Bspl(space)
    nanIdx = isnan(rBspl)
    rBspl[nanIdx] = 0
    return rBspl

#knotpoint construction (linearly spaced!)
def knotSeq(order,totalknots,RMIN,RMAX):
    btweenpts = totalknots-(2*(order-1))
    RMINARR = zeros(order-1)
    RMAXARR = RMAX*ones(order-1)
    RBTWN = linspace(0,RMAX,btweenpts) #The sequence between RMIN and RMAX
    t = append(RMINARR, RBTWN)
    t = append(t, RMAXARR)
    return t


