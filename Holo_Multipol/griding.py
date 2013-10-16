"""
   griding.py
   
   Using a predefined grid determine the volume element for every theta and every phi.
   
   -> In this program it is asumed that the grid is made up of different circles such that dtheta is constant and we can determine a dphi for every theta.
"""

__author__ = "Marion Baumgartner (marion.baumgartner@gmail.com)"
__date__ = "$Date: 14/10/2013 $"

from numpy import *

#Global variables
gridfile="oldinp.itp"
outputfile="dOmega.dat"
newFile="newinp.itp"


grid=list()


def makeGrid():
    """
    From the file oldinp.itp a grid created with a dphi for every point
    The dtheta that is calculated at the end is the same for all points
    """
    
    theta, phi=loadtxt(gridfile, usecols=(1,2), unpack=True)
    
    theta=list(theta)
    phi=list(phi)
    thetaTemp=list(set(theta))  #get rid of multiple values in the list so we have all different theta values    
    
    for i in range(len(theta)):
        grid.append([theta[i],phi[i],360/theta.count(theta[i]),0,0])
    dtheta=theta[0]-theta[1]
    if(dtheta<0):
        return -1*dtheta;
    else:
        return dtheta


def fitNewToOldGrid():
    """
    The function takes the nwe input data and fits it in to the grid that is given through the input file oldinp.itp
    
    The measured value of points that lay coles within a certain range is averaged and the value of the average is set the corresponding point of the given grid.
    In other words a two dimensional binning is done with the data given
    """
    gNew, thetaNew, phiNew=loadtxt(newFile, usecols=(0,1,2), unpack=True)
    gNew=list(gNew)
    thetaNew=list(thetaNew)
    phiNew=list(phiNew)
    
    if(len(grid)==0):
        dtheta=makeGrid()
    
    for i in range(len(gNew)):
        place=findPhi(phiNew[i],thetaNew[i],dtheta)
        if(place!=-1):
            grid[place][3]=grid[place][3]+gNew[i]
            grid[place][4]=grid[place][4]+1
        else:
            print "point (%f, %f) not found" %(phiNew[i],thetaNew[i])
            continue
    
    for i in range(len(grid)):
        if(grid[i][4]!=0):
            grid[i][3]=grid[i][3]/grid[i][4]
        else:
            continue
 
   
def findPhi(phi, theta, dtheta):
    """
    for a given phi and theta find the corresponding point in the grid and return the index list index of this point
    @param phi      double determining azimutal angle
    @param theta    double determining polar angle
    @return         int index of corresponing point on grid or -1 if point can not be fited on grid
    """
    if(len(grid)==0):
        print "no values in grid"
        return 0
    for i in range(len(grid)):
        if((theta<grid[i][0]+dtheta and theta<=grid[i][0]-dtheta) and (phi<grid[i][1]+grid[i][2] and phi>=grid[i][1]-grid[i][2])):
            return i
        else:
            continue
    return -1 #the walue could not be squeezed in to the new grid
    

dtheta=0

def main():
    dtheta=fitNewToOldGrid()
    with open(outputfile, 'w\n') as outFile:
        outFile.write("g(theta, phi), theta, phi, dphi (values are in degrees)\n")
        for i in range(len(grid)):
            value=str("%f %f %f %f" %(grid[i][3],grid[i][0],grid[i][1],grid[i][2]))
            #print value
            outFile.write(value+"\n")
    

if __name__ == '__main__':
    main()