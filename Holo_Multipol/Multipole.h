#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "Data.h"
#include <boost/math/special_functions/legendre.hpp>

class Multipole:public Data//:public Data //this class now has acces to all the pubilc and protected members of Data
                              //for acces to private members of Data make frined function/class
{
    public:
    /**
    *Class constructor
    *@param
    */
    Multipole(std::string fileName, int lmax, int isym);

    /**
    *Class constructor
    *@param
    */
    virtual ~Multipole();

    /**
    *   calculation of the real a_lm1 and imaginary a_lm2 coefficiants of the multipol expansion of the measured function g(theta,phi)
    *   The coefficients are calculated using the following formula:
    *               a_lm=int(Y_lm*(theta,phi) g(theta,phi)dOmega)
    *                   dOmega=sin(theta) dphi dtheta
    */
    void multpl(double dphi, double dtheta);

    void expans();

    static inline double deg_to_rad(double deg);
    static inline double rad_to_deg(double rad);

    private:
    const int LMAX; //maximum number of multipole expansion (amount of coefficients)
    const int ISYM; ////what is isym?

    private:

    std::vector<std::vector<double> > alm1; //real expansion coefficients
    std::vector<std::vector<double> > alm2; //imaginary expansion coefficients
//To add elements to the two dimensional vector use alm.push_back(row)
//where row is a vector of double type containing the data that belonges in to this row
    //std::vector<double> gCalc; take this vector from the base clas
};




#endif // MULTIPOLE_H
