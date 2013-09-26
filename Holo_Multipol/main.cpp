#include <iostream>
#include <cmath>

#include "Multipole.h"
#include "Data.h"


int main()
{
    Multipole *file=new Multipole("oldinp.itp",50,1,3.73);
    std::cout<<"created instance\n";
    file->readData();
    file->multpl();
    file->expans();

    //std::cout<<"apofact is "<<apofct(file->messg[0][1],file->messg[0][1]);

    return 0;
}








