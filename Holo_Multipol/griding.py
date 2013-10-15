"""
   griding.py
   
   Using a predefined grid determine the volume element for every theta and every phi.
   
   -> In this program it is asumed that the grid is made up of different circles such that dtheta is constant and we can determine a dphi for every theta.
"""

__author__ = "Marion Baumgartner (marion.baumgartner@gmail.com)"
__date__ = "$Date: 14/10/2013 $"

from numpy import *


gridfile="oldinp.itp"
outputfile="dOmega.dat"


theta, phi=loadtxt(gridfile, usecols=(1,2), unpack=True)
theta=list(theta)
phi=list(phi)
thetaTemp=list(set(theta))  #get rid of multiple values in the list so we have all different theta values


with open(outputfile, 'w\n') as outFile:
    outFile.write("#index, thetavalue, dphi (values are in degrees)\n")
    for i in thetaTemp:
        value=str("%d %.1f %.1f" %(theta.index(i), i, 360/theta.count(i)))
        #print value
        outFile.write(value+"\n")
    

