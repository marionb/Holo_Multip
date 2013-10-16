

grid=list()

def filGrid():
    for i in range(10):
        grid.append([i,i*10, i*100])
        
        
def changeOdd():
    if(len(grid)==0):
        filGrid()
    
    for i in range(len(grid)):
        if(grid[i][0]%2==0):
            grid[i][0]=0
            
        
filGrid()
print grid
changeOdd()
print grid