#include <iostream>
#include <cmath>

#include "Param.h"
#include "Data.h"


int main()
{
    Data *file=new Data("input.dat");
    std::cout<<"created instance\n";
    file->readData();

    //std::cout<<"apofact is "<<apofct(file->messg[0][1],file->messg[0][1]);

    return 0;
}








