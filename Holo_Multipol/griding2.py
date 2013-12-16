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
    def __init__(self, grid, expgrid=0, maxCoeff=50,):
        self.alm=list() #list containig the coefficients. Complex numbers
        self.grid=grid
        self.lmax=maxCoeff
        if(expgrid==0):
            pass
        else:
            self.expansionGrid=expgrid
            
    def degToRad(self, angle):
        return angle*math.pi/180.0
    
    def radToDeg(self, angle):
        return angle*180.0/math.pi
    
    """def multi(self, inorm=0, isym=1):
        bnorm=1
        for l in range(0,self.lmax+1,2):
            for m in range(l+1):
                rint=0
                iint=0
                
                for gr in self.grid:
                    gthph=gr[3]
                    polar=self.degToRad(gr[0])
                    azimuthal=self.degToRad(gr[1])
                    dphi=self.degToRad(gr[2])
                    norm=(2*l+1)**(0.5)
                    self.dtheta= 0.035
                    rint+=math.sin(polar)*self.dtheta*dphi*gthph*(sfunc.sph_harm(m,l,polar,azimuthal).real)*norm
                    iint-=math.sin(polar)*self.dtheta*dphi*gthph*(sfunc.sph_harm(m,l,polar,azimuthal).imag)
                
                if(l==0 and m==0 and rint>=0.001):
                    bnorm=rint
                #print "bnorm= ", bnorm
        
                self.alm1.append([l,m,rint/bnorm])
                self.alm2.append([l,m,iint/bnorm])"""
                
    def multi(self,inorm=0, isym=1):
        """
        calculates the coefficients A_lm for themulti pole expansion
        """
        bnorm=1.0
        for l in range(0,self.lmax+1,2):
            temp=list()
            for m in range(l+1):
                #Do the alculation
                
                summ=0
                print "-----------------"
                for val in self.grid:
                    temp2=sfunc.sph_harm(m,l,val[1],val[0])
                    summ+=temp2.conjugate()*val[5]*val[3]
                    
                if(l==0):
                    bnorm=summ
                    print"normalization factor= ",bnorm
                
                temp.append([l,m,summ,summ/bnorm])
                print temp[m]
            self.alm.append(temp)
            
    def intencity(self, theta, phi):
        """
        calculate the intencity g(theta,phi) at a given point theta, phi using
            g(theta,phi)=sum_{0<=l<=l_max}[A_l0*Y_l0(theta,phi))+2*Re(sum_{0<m<=l}A_lm*Y_lm(theta,phi))]
        
        @param theta    polar angle
        @param phi      azimutal angle
        """
        summ=0
        summ0=0
        for l in range(0,self.lmax+1,2):
            summ0+=(self.alm[l/2][0][3]*sfunc.sph_harm(0,l,phi,theta)).real
            for m in range(1,l+1):
                #print "-----------------------------------"
                summ+=(self.alm[l/2][m][3] * sfunc.sph_harm(m,l,phi,theta)).real
                
        return summ0+2*summ
            
    def expand(self):
        """
        calculate the intecities g(theta,phi) for all points on a given Sgrid
        """
        for term in self.grid:
            i=self.index(term)
            self.grid[i][2]=self.intencity(term[0],term[1])
    
    
def main():
    print "main function is not mentioned to be used!"
    pass  
    
if __name__ == '__main__':
    main()
    
    
