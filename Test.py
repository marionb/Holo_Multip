import random as rnd
import numpy as np

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
         
filGrid()
print grid
#changeOdd()
#print grid
a=list()
for i in range(10):
    numX=rnd.randint(10,40)
    numY=rnd.randint(10,40)

    numX=numX/10.0
    numY=numY/10.0
    a.append(numY)
    
    print numX, numY
    count = findPhi(numY,numX,1,1)
    print count
    print grid[count]
    

aArray=np.array(a)
print type(aArray)
print aArray

print np.where((1.5<aArray and aArray>2.5))


