#include "DataProcess.h"
#include <iostream>

DataProcess::DataProcess()
{

}
DataProcess::~DataProcess()
{

}

double DataProcess::deg_to_rad(double deg)
{
    return M_PI/180.0;
}

double DataProcess::apofct(double theta)
{
    //what is thetamax?is it just the maximum of all thetas?
    double width =2.0; //where does this number come from?? is it always 2?
    return 1.0/(std::exp((theta-thmax+3.0*width)/width) +1.0);
}


