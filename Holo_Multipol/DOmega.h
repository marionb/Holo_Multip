#ifndef DOMEGA_H
#define DOMEGA_H

#include <vector>
//#include <cmath>
//#include <assert.h>
#include <utility>  // std::pair, std::make_pair
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "Data.h"

//------------------------------------------------------------------------------------------//
template <class T1, class T2, class T3>
struct Triple
{    //Structure of three accociated elements of different types -> making life easier
    T1 first;
    T2 second;
    T3 third;
    /**
    * Constructor with given element
    */
    Triple(T1 e1, T2 e2, T3 e3):first(e1),second(e2),third(e3)
    {}

    /**
    * Default consstructor
    */
    Triple()
    {}
};
//------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------//
template<class T1, class T2, class T3, class T4>
struct Quadruple
{
    T1 first;
    T2 second;
    T3 third;
    T4 forth;
    /**
    * Constructor with given element
    */
    Quadruple(T1 e1, T2 e2, T3 e3, T4 e4):first(e1),second(e2),third(e3), forth(e4)
    {}

    /**
    * Default consstructor
    */
    Quadruple()
    {}
};
//------------------------------------------------------------------------------------------//

typedef std::vector< Triple<int, double, double> > Delta_Angle;
typedef std::vector< Quadruple<int, double, double, double, double> BaseGrid; //int->count; double->theta; double->phi; double->measured value g(theta,phi)
typedef std::pair<double,double> Omega;

class DOmega
{
    public:
    /**
    *Constructor
    *@param gridFiele   srting with name of the file containing the grid coordinates
    */
    DOmega(std::string gridFile);

    /**
    * for a given phi and theta find the maching DOmega value
    *@param phi     double value of azimutal angle  -> angle must be in radias
    *@param theta   double value of polar angle     -> angle must be in radias
    *@param return  functin return value consists of a pair of double values dphi and dtheta
    */
    Omega get_DOmega(double phi, double theta);

    void fitToGrid(Data element);

    private:

    void readFile();

    void set_dtheta();

    const char* gridFile;

    Delta_Angle dphi; //containing all thetas and coresponding dphis (it is assumed that for one teta we have one dphi)
    BaseGrid gridmes;
    double dtheta; //is assumed to be the same through out the whole grid

    std::ifstream grid;
};

#endif
