
from numpy import *
import sys
import griding2 as gr2

def main():
    
    print "Running", sys.argv[0]
    theta=0;
    phi=0;
    l=0;
    m=0;
    
    if len(sys.argv) == 7:#the name of the base grid and the name of the data fiel has to be given
        g=float(sys.argv[1])
        theta=float(sys.argv[2])
        phi=float(sys.argv[3])
        dOmega=float(sys.argv[4])
        l=float(sys.argv[5])
        m=float(sys.argv[6])
       
    else:
        print "to many arguments given program takes no argument or two arguments(name of fiel with grid, name of file containing data)"
        return -1
    
    
    print "intencity= ", gr2.intencityFree(g,theta,phi,dOmega, l ,m)
    
if __name__ == '__main__':
    main()