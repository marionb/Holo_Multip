import random as rnd
import numpy as np
from scipy.special import sph_harm
import math

grid=list()

def filGrid():
    for j in range(5):
        for i in range(5):
            grid.append([j,i])
        
        
def changeOdd():
    global grid
    if(len(grid)==0):
        filGrid()
    
    for i in range(len(grid)):
        if(grid[i][0]%2==0):
            grid[i][0]=0
            
            
            
def findPhi(phi, theta, deltheta, deltaphi):
    """
    for a given phi and theta find the corresponding point in the grid and return the index list index of this point
    @param phi      double determining azimutal angle
    @param theta    double determining polar angle
    @return         int index of corresponing point on grid or -1 if point can not be fited on grid
    """
    global grid
    if(len(grid)==0):
        print "no values in grid"
        return -1
    #print "foorloop printing variablers"
    for i in range(len(grid)):
        #print "theta range %f ,%f"%(grid[i][0]*1.0-deltheta/2.0, grid[i][0]*1.0+deltheta/2.0)
        #print "phi range %f ,%f"%(grid[i][1]*1.0-deltaphi/2.0, grid[i][1]*1.0-deltaphi/2.0)
        if(((grid[i][0]*1.0-deltheta/2.0)<=theta and theta<(grid[i][0]*1.0+deltheta/2.0)) and ((grid[i][1]*1.0-deltaphi/2.0)<=phi and phi<(grid[i][1]*1.0+deltaphi/2.0))):
            return i
        else:
            continue
    return -1 #the walue could not be squeezed in to the new grid
         


def spH():
    for l in range(4):
        m=0
           
        for theta in np.arange(0,2*math.pi,0.01):
            for phi in np.arange(0,math.pi,0.01):
                print l ,m, theta, phi, sph_harm(m,l,theta,phi)
                    



spH()
