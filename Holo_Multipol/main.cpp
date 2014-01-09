#include <iostream>
#include <cmath>
#include <complex>

#include "Multipole.h"
#include "Data.h"

int main()
{
    //double grid=0.09;
    Multipole *file=new Multipole("oldinp.itp",10,1,3.73);
    std::cout<<"created instance\n";
    file->readData();
    file->multpl();
    //file->printAlm();
    file->expans();

    //file->writeData("testOut.tex");
    //file->doyzimage(grid,0);

    //std::cout<<"apofact is "<<apofct(file->messg[0][1],file->messg[0][1]);

    delete file;
    return 0;
}








