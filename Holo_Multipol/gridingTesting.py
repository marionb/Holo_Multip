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
import sys
import griding2 as gr2



class Grid:
    #Global variables
    def __init__(self, gfile="oldinp.itp", nfile="newinp.itp"):
        self.gridfile=gfile   # file containing the grid information
        self.outputfile=""
        self.newFile=nfile    # file containing the measured data
               
        self.grid=list()                    # contais the grid and bining information theta,phi,dphi,g(theta,phi),count,dOmega_i
        self.theta=list()                   # contains all theta values at same index as corresponding phi value
        self.phi=list()                     # conains all phi walues
        self.thetaTemp=list()               # contains all different theta values only once. this list is creard from redundancy
        self.dtheta=0                       # the dtheta value for one data file of measurmet is assumed to be the same
        self.dOmega=list()                  # contains the resulting value of dOmega(theta,dphi,dtheta)=dphi*dtheta*sin(theta) -> using the dependancy dphi(theta). the values are in radian
    
    def makeGrid(self):
        """
        From the file oldinp.itp a grid created with a dphi for every point
        The dtheta that is calculated at the end is the same for all points
        """
        print"------------------------------------------------"
        
        print "clearing grid and writing new grid acording to file data"
        del self.grid[:]
        del self.theta[:]
        del self.phi[:]
        del self.thetaTemp[:]
        del self.dOmega[:]
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
        
        self.dtheta=abs(self.thetaTemp[0]-self.thetaTemp[1]) #dtheta is the same for all data points
        
        #for i in range(len(self.thetaTemp)):
        #    print self.thetaTemp[i], self.theta.count(self.thetaTemp[i])
        
        for i in range(len(self.theta)):
            
            dPhi=360.0/self.theta.count(self.theta[i])
            Omega=dPhi*self.dtheta*sin((self.theta[i])*pi/180)#dPhi*pi/180*self.dtheta*pi/180*sin((self.theta[i])*pi/180)
            self.dOmega.append(Omega)
            #print "dOmega= ",dPhi*pi/180*self.dtheta*pi/180*sin((self.theta[i])*pi/180)
            self.grid.append([self.theta[i],self.phi[i],dPhi,0,0, Omega])

        print "grid processing done."
        print"------------------------------------------------"    
    
    def calcGrid(self, phi_ref, theta_step, phi_range, maxTheta=90):
        """
        calculate the phi and theta walues used for the grid
        @param Phi_ref          Density of the phi steps for all theta values; Should be in correspondace with theta_step
        @param Theta_step       Size of the Theta steps (1, 2 oder 3 degree usw.)
        @param phi_range        The range for the phi values (0<phi_range<=360)
        """
        if(theta_step>maxTheta or theta_step<=0 or 0>phi_range or phi_range>360):
            try:
                raise SyntaxError('Wrong input')
            except SyntaxError:
                print 'Error found: '
                raise

        print "By using this function all vaues from a prewiously existing grid will be owerwritten!"
        print"------------------------------------------------"
        print"calculating grid using constant solid angles density;"
        
        del self.grid[:]
        del self.theta[:]
        del self.phi[:]
        del self.thetaTemp[:]
        del self.dOmega[:]
 
        self.dtheta=theta_step
        th=0
        #the first point in the grid is 0, We append this point here to make shure that it is in the grid
        self.grid.append([0,0,360,0,0,0])
        self.phi.append(0)
        self.theta.append(0)
        th+=theta_step
        while(th<=maxTheta):
            amount=int(phi_range*sin(th*pi/180)*phi_ref/theta_step) #function for constant solid angle density rounded to the neares whole number and converted int an integer
            #print th, amount
            if(amount!=0):
                dphi=360.0/amount
            else:
                dphi=0
            #print "amount, th dphi ",amount, th, dphi
            for i in range(int(amount)):
                #self.grid.append([th,i*dphi,dphi,0,0])
                self.theta.append(th)
                self.phi.append(i*dphi)
                self.dOmega.append(dphi*pi/180*self.dtheta*pi/180*sin((th)*pi/180)) #in radian
                self.grid.append([th,i*dphi,dphi,0,0,self.dOmega[i]])
                #print (th,i*dphi,dphi,0,0)
            th+=theta_step
            
        self.thetaTemp=list(set(self.theta))  #get rid of multiple values in the list so we have all different theta values
        gFile="grid.dat"
        print "calculated grid (used for the rest of the program) can be found in ",gFile
        with open(gFile, 'w\n') as outFile:
            for i in self.grid:
                value=str("%f %f" %(i[0],i[1]))
                outFile.write(value+"\n")
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
            print "Data fitted taken from", self.newFile
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
            #print "theta range %f, %f, %f"%(t-delYVal/2.0, t, t+delYVal/2.0)
            if(t-delYVal/2<=yVal and yVal<t+delYVal/2):
                findTheta=1
                index=self.theta.index(t)
                count=self.theta.count(t)
                #print "theta range %f, %f, %f"%(t-delYVal/2.0, t, t+delYVal/2.0)
                
                #calculate first interval to check
                phi1=self.grid[index][1]*1.0-self.grid[index][2]/2.0
                phi2=self.grid[index][1]*1.0+self.grid[index][2]/2.0
                
                for i in range(index,index+count):
                    phirange=[(phi1),(phi2)]
                    #print phirange
                    if((phirange[0]<=xVal and xVal<phirange[1]) or (phirange[1]>360 and (0<=xVal and xVal<phirange[1]%360)) or (phirange[0]<0 and (360+phirange[0]<=xVal and xVal>phirange[1]))):
                        return i
                    else:
                        #calculate new intervall
                        phi1=phi2
                        phi2=phi2+self.grid[i][2]
                        
            else:
                continue
        return -1
    
    
    def checkRange(self,ymin,y,ymax):
        """
        function checks if a given value lies within a cerain interval
        ==> IF y in [ymin,ymax[
            THEN return true
            ELSE return false
            
        @param ymin the minimum of the interval
        @param y    the value to be cheked
        @param ymax the maximum of the interval
        """
        if(ymin<=y and y<ymax):
            return True
        else:
            return False 
    
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
            #outFile.write("#g(theta, phi), theta, phi, dphi (values are in degrees)\n")
            for i in range(len(self.grid)):
                value=str("%f %f %f %f" %(self.grid[i][3],self.grid[i][0],self.grid[i][1],self.grid[i][5],))
                #print value
                outFile.write(value+"\n")
           # outFile.write("dtheta=%f"%self.dtheta)

    def OmegaArea(self):
        integral=0
        for i in self.dOmega:
            #print i
            integral+=i
        return integral
    
    def plotGrid(self):
        """
        The functin takes the points in given in the grid and plots them in a figure
        """
        import pylab as pyl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.axes3d import Axes3D
        tempx=[i*pi/180.0 for i in self.phi]
        tempy=[i*pi/180.0 for i in self.theta]
        print "(", self.phi[0],",",self.theta[0],")"
        print "(", tempx[0],",",tempy[0],")"
        x=cos(tempx)*sin(tempy)
        y=sin(tempx)*sin(tempy)
        z=cos(tempy)        
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c="r")
        plt.show()
        
        
    def writeGrid(self, gFile):
        print "calculated grid (used for the rest of the program) can be found in ",gFile
        with open(gFile, 'w\n') as outFile:
            for i in self.grid:
                value=str("%f %f %f %f %f %f" %(i[0],i[1],i[2],i[3],i[4],i[5]))
                outFile.write(value+"\n")
 
def main():
    
    print "Running", sys.argv[0]
    gfile=""
    nfile=""
    
    if len(sys.argv) == 3:#the name of the base grid and the name of the data fiel has to be given
        gfile=sys.argv[1]
        nfile=sys.argv[2]  
    elif len(sys.argv) == 1:
        gfile="oldinp.itp"
        nfile="newinp.itp"
    else:
        print "to many arguments given program takes no argument or two arguments(name of fiel with grid, name of file containing data)"
        return -1
    
    print "using grid information from  ",gfile
    print "using data from              ",nfile
    
    newGrid= Grid(gfile, nfile)
    newGrid.makeGrid()
#    newGrid.calcGrid(1, 2, 360)
    print "dOmega=", newGrid.OmegaArea()
    newGrid.fitNewToOldGrid()
    newGrid.writeGrid("temp.dat")
    #newGrid.plotGrid()
    
#    calc = gr2.Calc(newGrid.grid,0)
    
#    calc.multi()
#    calc.expand()
#    calc.writeData("calc.dat")
    
    
    #print('\a')
    #sys.stdout.write('\a')
    #sys.stdout.flush()
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
    
    
