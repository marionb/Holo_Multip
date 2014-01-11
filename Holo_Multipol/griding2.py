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
import scipy.special as sfunc
import math

class Calc:
    #Global variables
    def __init__(self, grid, maxCoeff=50,):
        self.alm=list() #list containig the coefficients. Complex numbers
        self.grid=grid  #needs to be of format theta phi g(theta,phi)
        #self.endGrid    #needs to be of format theta phi. g_calc(theta,phi) will be apended
        self.lmax=maxCoeff
            
    def degToRad(self, angle):
        return angle*math.pi/180.0
    
    def radToDeg(self, angle):
        return angle*180.0/math.pi
                
    def multi(self,inorm=0, isym=1):
        """
        calculates the coefficients A_lm for themulti pole expansion
        """
        print "using l_max=",self.lmax
        bnorm=1.0
        for l in range(0,self.lmax+1,2):
            temp=list()
            for m in range(l+1):
                #Do the calculation
                
                summ=0
                #print "-----------------"
                for val in self.grid:
                    temp2=sfunc.sph_harm(m,l,val[1],val[0])
                    summ+=temp2.conjugate()*val[3]*val[5] #Y*_lm(theta,phi)*g(thea,phi)*sin(theta)dphi*dtheta
                    
                if(l==0 and m==0):
                    bnorm=summ
                    print"normalization factor= ",bnorm
                
                temp.append([l,m,summ,summ/bnorm])
                #print temp[m]
            self.alm.append(temp)
        self.writeCoeff("PYcoeff.dat")
        print "calcualted coefficients"
        print "-------------------------------------------"
        
    def intencity(self, theta, phi):
        """
        calculate the intencity g(theta,phi) at a given point theta, phi using
            g(theta,phi)=sum_{0<=l<=l_max}[A_l0*Y_l0(theta,phi))+2*Re(sum_{0<m<=l}A_lm*Y_lm(theta,phi))]
        
        @param theta    polar angle
        @param phi      azimutal angle
        """
        summ0=0
        summ=complex(0,0)
        for l in range(0,self.lmax+1,2):
            summ0+=(self.alm[l/2][0][3]*sfunc.sph_harm(0,l,phi,theta)).real
            #print self.alm[l/2]
            for m in range(1,l+1):
                #print "-----------------------------------"
                summ+=self.alm[l/2][m][3]* sfunc.sph_harm(m,l,phi,theta) 
        return summ0+2*summ.real
            
    def expand(self):
        """
        calculate the intecities g(theta,phi) for all points on a given Sgrid
        """
        print "expanding new function using grid ..."
        for term in self.grid:
            i=self.grid.index(term)
            self.grid[i][2]=self.intencity(term[0],term[1])
            
        print "function has been expanded"
        print "--------------------------------------"
        
    def writeData(self, dataFiel):
        """
        write information from grit in to dataFile
        """
        print "Calculated Data can be found in ",dataFiel
        with open(dataFiel, 'w\n') as outFile:
            for i in self.grid:
                value=str("%f %f %f" %(i[0],i[1],i[2]))
                outFile.write(value+"\n")
    
    def writeCoeff(self, dataFiel):
        """
        write information from Coeff in to dataFile
        """
        print "Calculated Coefficients can be found in ",dataFiel
        with open(dataFiel, 'w\n') as outFile:
            value1=str(self.lmax)
            outFile.write(value1+"\n")
            
            for i in self.alm:
                for h in i:
                    value=str("%f %f %f %f" %(h[0], h[1], h[3].real, h[3].imag))
                    outFile.write(value+"\n")            


def intencityFree(g,theta, phi, dOmega, l ,m):
        """
        calculate the intencity g(theta,phi) at a given point theta, phi using
            g(theta,phi)=sum_{0<=l<=l_max}[A_l0*Y_l0(theta,phi))+2*Re(sum_{0<m<=l}A_lm*Y_lm(theta,phi))]
        
        @param theta    polar angle
        @param phi      azimutal angle
        """
        temp2=sfunc.sph_harm(m,l,phi,theta)   
        return temp2.conjugate()*g*dOmega


def main():
    print "main function is not ment to be used!"
    pass  
    
if __name__ == '__main__':
    main()
    
    
