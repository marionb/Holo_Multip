#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <assert.h>
#include "Data.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/bessel.hpp>

class Multipole:public Data //:public Data //this class now has acces to all the pubilc and protected members of Data for acces to private members of Data make frined function/class
{
    public:
    /**
    *Class constructor
    *@param
    */
    Multipole(std::string fileName, int lmax, int isym, double ekin);

    /**
    *Class constructor
    *@param
    */
    virtual ~Multipole();

    /**
    *returne the member variable radimg
    */
    std::vector<double> get_radimg();

    /**
    *   calculation of the real a_lm1 and imaginary a_lm2 coefficiants of the multipol expansion of the measured function g(theta,phi)
    *   The coefficients are calculated using the following formula:
    *               a_lm=int(Y_lm*(theta,phi) g(theta,phi)dOmega)
    *                   dOmega=sin(theta) dphi dtheta
    */
    void multpl();

    void expans();

    /**
    *calculateion of a radial image curve of the electron wave field near the photoemitter.
    *Equation (3) from A.Stucke et al. (1992) is used.
    *@param alpha, beta     double precission; define the radial curve; the angles are given within the same polar coordinates frame as the data.
    *The radial grid is in angstrom with 0.1\AA spacing
    *Function is adopted from the fortran Program by Juerg Osterwalde writen in 1993
    */
    void holorad(double, double);

    /**
    *Calculates a two dimensional image of the electron wave field near the photoemitter.
    *Equation (3) from A.Stucke et al. (1992) is used.
    *Ther radial grid is in angstroem with grid spacing in Angstroe
    *Function is adopted from the fortran Program by Juerg Osterwalde writen in 1993
    */
    void doyzimage(double grid, int xyz);

    void printAlm();

    private:
    const int LMAX; //maximum number of multipole expansion (amount of coefficients)
    const int ISYM; ////what is isym?

    double k; //2*pi*sqrt(ekin/150)

    std::vector<std::vector<double> > alm1; //real expansion coefficients
    std::vector<std::vector<double> > alm2; //imaginary expansion coefficients
    //To add elements to the two dimensional vector use alm.push_back(row)
    //where row is a vector of double type containing the data that belonges in to this row
    //std::vector<double> gCalc; take this vector from the base clas

    /**
    *calculate (-1)^exp with exp in [1,2,3,4,...]
    *@param exp is a positiv or negative integer
    */
    int vorz(int);


    inline double innerSumm(int l, double alpha, double beta);

    /**
    *calculation of the proper polar angle within the yz-plane.
    *All polar angles are positive; negative y values should be accounted for by adding PI to the value of phi.
    *@return polar angle in radian (double precission)   ->  if the function returnes -1 there calculation was corrupded or the function is wrong
    */
    double calcth(double y, double z);

    /**
    *calculation of the proper azimutal angle within the yz-plane.
    *@return azimutal angle in radian (double precission)   ->  if the function returnes -1 there calculation was corrupded or the function is wrong
    */
    double calcphi(double x, double y);
};


#endif // MULTIPOLE_H
