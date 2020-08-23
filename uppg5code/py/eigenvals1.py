from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as LA #linear algebra
from bspline_knot_functs import B,knotSeq #Bspline and knotsequence construction
from Hmatrix import Hmat,H_GQint
from Bmatrix import Bmat,B_GQint
from globals import * #imports l,Z,N,k,rmin,rmax - the global parameters.

Hh = Hmat(N,k,rmin,rmax,l,Z)
Bb = Bmat(N,k,rmin,rmax)
print diag(Bb)
print Bb
#print Hh[-1][-1], Bb[-1][-1] #SEEMS TO WORK
