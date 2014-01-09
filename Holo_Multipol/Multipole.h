#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <assert.h>
#include <complex>

#include "Data.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/bessel.hpp>

class Multipole:public Data //:public Data //this class now has acces to all the pubilc and protected members of Data for acces to private members of Data make frined function/class
{
    public:
    /**
    *Class constructor
    *@param fileName    name of the input file
    *@param lmax        TODO Describe
    *@param isym        TODO Describe
    *@param ekin        kinetic enery
    */
    Multipole(std::string fileName, int lmax, int isym, double ekin);

    /**
    *getter method returning the number of expansion parameters
    */
    const int getLMAX();

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
    void multpl();


    /**
    *Calcualation of the angular distribution g(theta,phi) from multipole expansion coefficients
    *alm1 and alm2 (real and imaginary parts of alm). The equation (3.58) from Jackson, p. 100 is used.
    *g is real so that we can use only positive m values and the equation becomes
    *
    *g = sum_l ( a_l0 * Y_l0 + 2 * Re sum_m (a_lm * Y_lm) ).
    *
    *The sum_l runs from 0 to lmax and sum_m from 1 to l.
    *Adapted from the corresponding Fortran subroutin by Juerg Osterwalder, universite de fribourg, 2.7.93
    */
    void expans();

    /**
    *calculateion of a radial image curve of the electron wave field near the photoemitter.
    *Equation (3) from A.Stucke et al. (1992) is used.
    *@param alpha, beta     double precission; define the radial curve; the angles are given within the same polar coordinates frame as the data.
    *The radial grid is in angstrom with 0.1\AA spacing
    *Function is adopted from the fortran Program by Juerg Osterwalde writen in 1993
    */
    //void holorad(double, double);

    /**
    *Calculates a two dimensional image of the electron wave field near the photoemitter.
    *Equation (3) from A.Stucke et al. (1992) is used.
    *Ther radial grid is in angstroem with grid spacing in Angstroe
    *Function is adopted from the fortran Program by Juerg Osterwalde writen in 1993
    */
    //void doyzimage(double grid, int xyz);

    /**
    *Calculates a two-dimensional image of the electron wave field near the photoemitter.
    *Equation (3) from A.Stucke et al. (1992) is used.
    *the radial grid is in angstroem with grid angstroem spacing.
    *Function is adapted from the fortran program by Juerg Osterwalder written in 6.7.93
    */
    //void doxyzimage(double grid);

    /**
    *Takes a 2D image array and maps it into a 2d image array of integer numbers ranging from 0 to 255.
    *Function makes input in to 'image program more convenient for mac.
    *Calculates r * |img(xy)|**2 to produce images.
    *adapted from fortran program by Juerg Osterwalder, universite de fribourg, 12.7.93
    */
    //void scaleimage(double grid);

    /**
    *Function performes a gausian smoothing of the image function
    */
    //void smooth(double grid);

    /**
    *print the real coefficients that were calculated within this program
    */
    //void printAlm();

//=======================================================================================//
    private:
    const int LMAX; //maximum number of multipole expansion (amount of coefficients)
    const int ISYM; ////what is isym?

    const double k; //2*pi*sqrt(ekin/150)

    std::vector<std::vector<std::complex<double> > > alm; //expansion coefficients

    //To add elements to the two dimensional vector use alm.push_back(row)
    //where row is a vector of double type containing the data that belonges in to this row
    //std::vector<double> gCalc; take this vector from the base clas

    /**
    *calculate (-1)^exp with exp in [1,2,3,4,...]
    *@param exp is a positiv or negative integer
    */
    //int vorz(int);


    /**
    *calculate the intencity g(theta,phi) at a given point theta, phi using
    *    g(theta,phi)=sum_{0<=l<=l_max}[A_l0*Y_l0(theta,phi))+2*Re(sum_{0<m<=l}A_lm*Y_lm(theta,phi))]

    *@param theta    polar angle
    *@param phi      azimutal angle
    */
    inline double intencity(double theta, double phi);

    /**
    *calculation of the proper polar angle within the yz-plane.
    *All polar angles are positive; negative y values should be accounted for by adding PI to the value of phi.
    *@return polar angle in radian (double precission)   ->  if the function returnes -1 there calculation was corrupded or the function is wrong
    */
    //double calcth(double y, double z);

    /**
    *calculation of the proper azimutal angle within the yz-plane.
    *@return azimutal angle in radian (double precission)   ->  if the function returnes -1 there calculation was corrupded or the function is wrong
    */
    //double calcphi(double x, double y);

    static const int MAX_COEFF=100;//100 is the maximum amunt of coefficients that can be calculated
};


#endif // MULTIPOLE_H
