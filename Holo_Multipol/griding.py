"""
griding.py

The script takes a given two dimensional grid of points and a set of Data points.
The Data points are mached in to the grid and the measured values at every point is averaged over the number of points that correspont to one grid point.
(In other words a two dimensional bining is doen using a given grid.)
The Grid that is used is a Data set of a circular measurement with the original Detector.
     -> In this program it is assumed that the grid is made up of different circles such that dtheta is constant and we can determine a dphi for every theta.

TODO:Make it moer user friendly: EG place option to choos the name of input files -> use sys...
"""

__author__ = "Marion Baumgartner (marion.baumgartner@gmail.com)"
__date__ = "$Date: 14/10/2013 $"

from numpy import *

class Grid:
    #Global variables
    def __init__(self):
        self.gridfile="oldinp.itp"   # file containing the grid information
        self.outputfile=""
        self.newFile="newinp.itp"    # file containing the measured data
               
        self.grid=list()                    # contais the grid and bining information theta,phi,dphi,g(theta,phi),count
        self.theta=list()                   # contains all theta values at same index as corresponding phi value
        self.phi=list()                     # conains all phi walues
        self.thetaTemp=list()               # contains all different theta values oned
        self.dtheta=0
    
    def makeGrid(self):
        """
        From the file oldinp.itp a grid created with a dphi for every point
        The dtheta that is calculated at the end is the same for all points
        """
        print"------------------------------------------------"
        try:
            self.theta, self.phi=loadtxt(self.gridfile, usecols=(1,2), unpack=True)
            print "Beginning grid processing"
        except:
            print "data file could not be opend. Check that ",gridfile," exists in runing directory and run again."
            print"------------------ERROR-------------------------"
            print "will exit program!"
            print"------------------------------------------------"
            
        self.theta=list(self.theta)
        self.phi=list(self.phi)
        self.thetaTemp=list(set(self.theta))  #get rid of multiple values in the list so we have all different theta values    
        
        for i in range(len(self.theta)):
            self.grid.append([self.theta[i],self.phi[i],360.0/self.theta.count(self.theta[i]),0,0])
            
        self.dtheta=abs(self.thetaTemp[0]-self.thetaTemp[1]) #dtheta is the same for all data points

        print "grid processing done."
        print"------------------------------------------------"
    
    
    def fitNewToOldGrid(self):
        """
        The function takes the nwe input data and fits it in to the grid that is given through the input file oldinp.itp
        
        The measured value of points that lay coles within a certain range is averaged and the value of the average is set the corresponding point of the given grid.
        In other words a two dimensional binning is done with the data given
        
        global grid
        global newFile
        global theta
        global phi"""
        
        if(len(self.grid)==0 or self.dtheta==0):
            self.makeGrid()
            
        print"------------------------------------------------"
        try:
            gNew, thetaNew, phiNew=loadtxt(self.newFile, usecols=(0,1,2), unpack=True)
            print "Beginning data processing; this will take a while"       
            
        except:
            print "data file could not be opend. Check that ",newFile," exists in runing directory and run again."
            print"------------------ERROR-------------------------"
            print "will exit program!"
            print"------------------------------------------------"
        gNew=list(gNew)
        thetaNew=list(thetaNew)
        phiNew=list(phiNew)
        
        
        #    H,thetaEdge,phiedge=histogram2d(thetaNew,phiNew,bins=[theta,phi],weights=gNew)
        
        for i in range(len(gNew)):
            #print"--------------------"
            #print phiNew[i], thetaNew[i]
            place=self.findPhi(phiNew[i],thetaNew[i])
            if(place!=-1):
                self.grid[place][3]=self.grid[place][3]+gNew[i]
                self.grid[place][4]=self.grid[place][4]+1
            else:
                print "point %d (%f, %f) not found" %(i, phiNew[i],thetaNew[i])
                continue
        
        for i in range(len(self.grid)):
            if(self.grid[i][4]!=0):
                self.grid[i][3]=self.grid[i][3]/self.grid[i][4]
            else:
                continue
        print "mached new data in to base grid."
        print"------------------------------------------------"
    
       
    def findPhi(self, xVal, yVal):
        """
        for a given phi and theta find the corresponding point in the grid and return the index list index of this point
        @param xVal      double determining azimutal angle
        @param yVal      double determining polar angle
        @return          int index of corresponing point on grid or -1 if point can not be fited on grid
        global grid
        global theta"""
         
        
        if(len(self.grid)==0):
            print "no values in grid"
            return -1
        
        delYVal=self.dtheta
        
        for t in self.thetaTemp:
            #print "theta range %f, %f, %f"%(t*1.0-delYVal/2.0, t, t*1.0+delYVal/2.0)
            if(t-delYVal<=yVal and yVal<t+delYVal):
                index=self.theta.index(t)
                count=self.theta.count(t)
                return index
                
                for i in range(index,count):
                    phirange=[(selg.grid[i][1]*1.0-self.grid[i][2]/2.0),(self.grid[i][1]*1.0+self.grid[i][2]/2.0)]
                    if(self.grid[i][1]*1.0-self.grid[i][2]/2.0<=xVal and xVal<self.grid[i][1]*1.0+self.grid[i][2]/2.0):
                        return i
                    else:
                        continue
                        
            else:
                continue
        return -1
        
    
    def fortranOut(self):
        self.outputfile="FortranOut.dat" # the output Data file
        with open(self.outputfile, 'w\n') as outFile:
            for i in range(len(self.grid)):
                value=str("%f %f %f" %(self.grid[i][3],self.grid[i][0],self.grid[i][1]))
                #print value
                outFile.write(value+"\n")
    
    def cPPOut(self):
        
        self.outputfile="CPPOut.dat" # the output Data file
        with open(self.outputfile, 'w\n') as outFile:
            outFile.write("#g(theta, phi), theta, phi, dphi (values are in degrees)\n")
            for i in range(len(self.grid)):
                value=str("%f %f %f %f" %(self.grid[i][3],self.grid[i][0],self.grid[i][1],self.grid[i][2]))
                #print value
                outFile.write(value+"\n")
            outFile.write("dtheta=%f"%self.dtheta)
    
def main():
    print "Runing griding.py"
    newGrid= Grid()
    newGrid.fitNewToOldGrid()
    s=-1
    while True:
        s = int(raw_input('Type 0 for Fortan output; 1 for c++ output '))
        if(s==0):
            newGrid.fortranOut()
            print"----------------------------------------------------"
            print "writing output for Fortran code to ",newGrid.outputfile
            print"----------------------------------------------------"
            return
        if(s==1):
            newGrid.cPPOut()
            print"----------------------------------------------------"
            print "writing output for C++ code to ",newGrid.outputfile
            print"----------------------------------------------------"
            return
        
    
        
    
if __name__ == '__main__':
    main()
    
    