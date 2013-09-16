#ifndef DATAPROCESS_H
#define DATAPROCESS_H

#define _USE_MATH_DEFINES
#include <cmath>

class DataProcess:public Data //this class now has acces to all the pubilc and protected members of Data
                              //for acces to private members of Data make frined function/class
{
    public:
    /**
    *Class constructor
    *@param
    */
    DataProcess();

    /**
    *Class constructor
    *@param
    */
    virtual ~DataProcess();

    //void multpl();

    //void expans();

    //double calcth(double y, double z);

    /**
    *calculation of a fermi function like apodizing function
    *PRE:   double theta and double thmax -> values for which te function is calculated
    *POST:  double with the value of the apodizing function
    */
    double apofct(double theta, double thmax);

    static double deg_to_rad(double deg);

    private:

    protected:


};




#endif // DATAPROCESS_H
