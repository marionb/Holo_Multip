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
typedef std::vector<std::vector<std::complex<dataType> > > AlmType;

class Multipole:public Data //this class now has acces to all the pubilc and protected members of Data for acces to private members of Data make frined function/class
{
    public:
    /**
    *Class constructor
    *@param fileName    name of the input file
    *@param apo         determine wther the apodization function is caled on the walues red from the data fiel
    *@param inorm       fag to determine wether all alm's are used. If this is grater than 0 the contribution of all alm with l<= inorm are cut/set to zero
    *@param lmax        TODO Describe
    *@param isym        TODO Describe
    */
    Multipole(std::string fileName, bool apo, int inorm, int lmax, int isym);

    /**
    *getter method returning the number of expansion parameters
    */
    const int getLMAX();

    /**
    *Class destructor
    *@param -
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
    *print the real coefficients that were calculated within this program
    */
    void writeAlm(std::string almFile);

    /**
    * read alm values from a data file
    *@param almFile     file name containig the alm with the folowing structure
    *       #l  m   A_lm
    *        0  0   (1,0)
    *        2  0   (0.38151,0)
    *        2  1   (-0.0132918,-0.0278365)
    *        2  2   (0.00317809,-0.000486047)
    *               ...     ...
    * a line containing a # is considerd to contain comment and is not read
    * the Complex number format of the alm are given as (real, complex)
    */
    void readAlm(std::string almFile);





//=======================================================================================//
    private:
    static const int MAX_COEFF=100;//100 is the maximum amunt of coefficients that can be calculated
    const int LMAX; //maximum number of multipole expansion (amount of coefficients)
    const int ISYM;
    const int INORM;

    //const dataType k; //2*pi*sqrt(ekin/150)

    AlmType alm; //expansion coefficients

    //To add elements to the two dimensional vector use alm.push_back(row)
    //where row is a vector of dataType type containing the data that belonges in to this row
    //std::vector<dataType> gCalc; take this vector from the base clas

    /**
    * The function uses the flag INORM and sets all alm with l<=INORM to zero
    * This is done to remove alm that have a large oszillation which influences the calculated g(theta,phi) function in an unwanted way
    */
    void removeAlmBackground();


    /**
    *calculate the intencity g(theta,phi) at a given point theta, phi using
    *    g(theta,phi)=sum_{0<=l<=l_max}[A_l0*Y_l0(theta,phi))+2*Re(sum_{0<m<=l}A_lm*Y_lm(theta,phi))]

    *@param theta    polar angle
    *@param phi      azimutal angle
    */
    dataType intencity(dataType theta, dataType phi);



};


#endif // MULTIPOLE_H
