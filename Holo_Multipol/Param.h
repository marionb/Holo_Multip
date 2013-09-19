#ifndef PARAM_H
#define PARAM_H

#include <fstream>
#include <iostream>
#include <limits>

class Param
{
    public:
        Param();
        virtual ~Param();
        void readParafiel();

    private:
        std::ifstream param;
        char parfile [];
        int lmax; //max amount of multiole expansions
        double ekin; //electron kinetic energy
        std::string restParams [];
};

#endif // PARAM_H
