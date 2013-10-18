"""
   griding.py
   
   The script takes a given two dimensional grid of points and a set of Data points.
   The Data points are mached in to the grid and the measured values at every point is averaged over the number of points that correspont to one grid point.
   (In other words a two dimensional bining is doen using a given grid.)
   The Grid that is used is a Data set of a circular measurement with the original Detector.
        -> In this program it is assumed that the grid is made up of different circles such that dtheta is constant and we can determine a dphi for every theta.
"""

__author__ = "Marion Baumgartner (marion.baumgartner@gmail.com)"
__date__ = "$Date: 14/10/2013 $"

from numpy import *

#Global variables
gridfile="oldinp.itp"   # file containing the grid information
outputfile=""
newFile="newinp.itp"    # file containing the measured data


grid=list()             # contais the grid and bining information theta,phi,dphi,g(theta,phi),count


def makeGrid():
    """
    From the file oldinp.itp a grid created with a dphi for every point
    The dtheta that is calculated at the end is the same for all points
    """
    global grid
    global gridfile
    theta, phi=loadtxt(gridfile, usecols=(1,2), unpack=True)
    
    theta=list(theta)
    phi=list(phi)
    thetaTemp=list(set(theta))  #get rid of multiple values in the list so we have all different theta values    
    
    for i in range(len(theta)):
        grid.append([theta[i],phi[i],360/theta.count(theta[i]),0,0])
        
    dtheta=theta[0]-theta[theta.count(theta[0])+1] #dtheta is the same for all data points
    
    if(dtheta<0):
        return -1*dtheta;
    else:
        return dtheta


def fitNewToOldGrid(dpolar):
    """
    The function takes the nwe input data and fits it in to the grid that is given through the input file oldinp.itp
    
    The measured value of points that lay coles within a certain range is averaged and the value of the average is set the corresponding point of the given grid.
    In other words a two dimensional binning is done with the data given
    """
    global grid
    global newFile
    gNew, thetaNew, phiNew=loadtxt(newFile, usecols=(0,1,2), unpack=True)
    gNew=list(gNew)
    thetaNew=list(thetaNew)
    phiNew=list(phiNew)
    
    if(len(grid)==0 or dpolar==0):
        dpolar=makeGrid()
        
    for i in range(len(gNew)):
        #print"--------------------"
        #print phiNew[i], thetaNew[i]
        place=findPhi(phiNew[i],thetaNew[i],dpolar)
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
    return dpolar
 
   
def findPhi(phi, theta, deltheta):
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
    
    for i in range(len(grid)):
        #print "theta range %f ,%f"%(grid[i][0]*1.0-deltheta/2.0, grid[i][0]*1.0+deltheta/2.0)
        #print "phi range %f ,%f"%(grid[i][1]*1.0-deltaphi/2.0, grid[i][1]*1.0+deltaphi/2.0)
        if(((grid[i][0]*1.0-deltheta/2.0)<=theta and theta<(grid[i][0]*1.0+deltheta/2.0)) and ((grid[i][1]*1.0-grid[i][2]/2.0)<=phi and phi<(grid[i][1]*1.0+grid[i][2]/2.0))):
            return i
        else:
            continue
    return -1 #the walue could not be squeezed in to the new grid
    

def fortranOut(grid):
    global outputfile
    outputfile="FortranOut.dat" # the output Data file
    with open(outputfile, 'w\n') as outFile:
        for i in range(len(grid)):
            value=str("%f %f %f" %(grid[i][3],grid[i][0],grid[i][1]))
            #print value
            outFile.write(value+"\n")

def cPPOut(grid, dtheta):
    global outputfile
    outputfile="CPPOut.dat" # the output Data file
    with open(outputfile, 'w\n') as outFile:
        outFile.write("#g(theta, phi), theta, phi, dphi (values are in degrees)\n")
        for i in range(len(grid)):
            value=str("%f %f %f %f" %(grid[i][3],grid[i][0],grid[i][1],grid[i][2]))
            #print value
            outFile.write(value+"\n")
        outFile.write("dtheta=%f"%dtheta)

def main():
    dtheta=0
    dtheta=fitNewToOldGrid(dtheta)
    s=-1
    while True:
        s = int(raw_input('Type 0 for Fortan output; 1 for c++ output '))
        if(s==0):
            fortranOut(grid)
            print "writing output for Fortran code to ",outputfile
            return
        if(s==1):
            cPPOut(grid, dtheta)
            print "writing output for C++ code to ",outputfile
            return
    

    

if __name__ == '__main__':
    main()