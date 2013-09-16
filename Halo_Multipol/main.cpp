#include <iostream>
#include <cmath>

#include "Param.h"
#include "Data.h"

int main()
{
    Data *file=new Data("input.dat");
    std::cout<<"created instance\n";
    file->readData();

    return 0;
}




/**
*calculation of a fermi function like apodizing function
*PRE:   double theta and double thmax -> values for which te function is calculated
*POST:  double with the value of the apodizing function
*/
/*double apofct(double theta, double thmax)
{
    double width =2.0; //wher does this number come from?? is it always 2?
	return 1.0/(std::exp((theta-thmax+3.0*width)/width) +1.0);
}*/



