#ifndef DATA_H
#define DATA_H

#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <vector>


class Data
{
    public:
        /**
        *Class constructor
        *@param input   char array of data file name
        */
        Data(std::string fileName);

        virtual ~Data();

        /**
        *Read data file and wirite the data in to the 2D array messg
        *The angles that are read are tranformed from degree to radian making the rest of the calculations easier
        *Also the largest theta angle is determined bythe function
        */
        void readData();

        /**
        * determine the amount of lines in the file
        * @param fileGiven  name of input file
        */
        static int getLineNum(std::string fileGiven);


        /**
        * write data from messgFinal in to an outpu file
        * @param output     name of the file to write data in
        */
        void writeData(std::string output);

        //------------------------------------------------------//
        //Static functions

        /**
        *convert an angle in degre to rad
        *@param deg double, value of the angle in degrees
        *@return    value of the angle in radians using double precission
        */
        static inline double deg_to_rad(double deg);

        /**
        *convert an angle in rad to degrees
        *@param rad double, value of the angle in radian
        *@return    value of the angle in degrees using double precission
        */
        static inline double rad_to_deg(double rad);

        //friend class Multipole;


    protected:

        /**
        *calculation of a fermi function like apodizing function
        *PRE:   double theta and double thmax -> values for which te function is calculated
        *POST:  double with the value of the apodizing function
        */
        inline void apofct();

        /**
        * calculates the chi function from the measured hologram
        * according to the procedure of chuck fadley, i.e. by
        * normalizing the data on each phi-circle to zero average.
        */
        void calcchi();


        const int MAXANGLES; //number of data points
        const int ROWNUMBER; //number of rows in the file read


//        bool iapo;
        double thmax;

//    protected:
        std::vector<std::vector<double> > messg;
        std::vector<std::vector<double> > calc;
        //double (*messg)[3]; //dynamic array containing input data -> this data needs to be accesible
        //double (*messgFinal)[3]; //dynamic array containing output data

    private:
        const char* infile; //name of file with input data
        std::ifstream dataF; //for input
        std::ofstream dataW; //for output
};

#endif // PARAM_H
