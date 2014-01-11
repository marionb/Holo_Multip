import scipy as sp
import numpy as np
import scipy.special as sfunc



def y_10(theta,phi):
    return complex(1,0)*0.5*(3/np.pi)**(0.5)* sp.cos(theta)

def y_11(theta,phi):
    return np.exp(complex(0,phi))*(-1/2.0)*(3/(2*np.pi))**0.5*sp.sin(theta)
def y_21(theta,phi):
    pass

def test():
    print 0.5*(1/np.pi)**0.5, sfunc.sph_harm(0,0,0,0)
    for theta in sp.arange(0,np.pi,np.pi/4.0):
        for phi in sp.arange(0,2*np.pi,np.pi/2.0):                
            print theta, phi,"===>", y_10(theta,phi)/sfunc.sph_harm(0,1,phi,theta),"    ", y_11(theta,phi)/sfunc.sph_harm(1,1,phi,theta)
            
test()