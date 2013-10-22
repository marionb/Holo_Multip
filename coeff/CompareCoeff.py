
import numpy as np

fortran="fortran.alm"
cpp="temp"

fReC,fImC=np.loadtxt(fortran, usecols=(2,3), unpack=True)
cReC,cImC=np.loadtxt(cpp, usecols=(2,3), unpack=True)


for i,j in zip(fReC, cReC):
    print abs(i-j), i, j
    