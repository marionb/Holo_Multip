#ifndef DATA_H
#define DATA_H

#include <fstream>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <vector>
#include "specialVector.h"

typedef float dataType;
typedef std::vector<dataType> oneVector;
typedef std::vector<std::vector<dataType> > twoVector; //a two dimensional vectors
typedef std::vector<SpecialVector<dataType> > twoVectorxNyZ; //a two dimensional vector with negative and positive y positions
typedef std::vector<std::vector<std::vector<dataType> > > threeVector;
typedef std::vector<SpecialVector<SpecialVector<dataType> > > threeVectorxNyZzZ;//a three dimensional vector with negative and positive y and z positions
class Data
{
    public:

        Data(std::string fileName, bool apo, bool inorm);
        /**
        *Class constructor
        *@param input   char array of data file name
        */
        virtual ~Data();

        /**
        *Read data file and wirite the data in to the 2D array messg
        *The angles that are read are tranformed from degree to radian making the rest of the calculations easier
        *Also the largest theta angle is determined bythe function
        */
        void readData();

        /**
        *Caling this function clears the vector grid and fills it with informaton taken from the file with the name given to the function
        *@param gridFile    string containing the file name with the grid
        *--> the file describen above has to contain two rows of float/double values
        *       first row       seccond row
        *       theta           phi
        */
        void readGrid(std::string);


        /**
        * write data from messgFinal in to an outpu file
        * @param output     name of the file to write data in
        */
        void writeData(std::string output);

        /**
        *Print the grid used for the calculation of g(theta,phi) on to the CMD-line
        */
        void printGrid();

        //------------------------------------------------------//
        //Static functions

        /**
        * determine the amount of lines in the file
        * @param fileGiven  name of input file
        */
        static int getLineNum(std::string fileGiven);

        /**
        *convert an angle in degre to rad
        *@param deg dataType, value of the angle in degrees
        *@return    value of the angle in radians using dataType precission
        */
        static inline dataType deg_to_rad(dataType deg);

        /**
        *convert an angle in rad to degrees
        *@param rad dataType, value of the angle in radian
        *@return    value of the angle in degrees using dataType precission
        */
        static inline dataType rad_to_deg(dataType rad);




    protected:

        const int MAXANGLES; //number of data points
        const int ROWNUMBER; //number of rows in the file read

        dataType thmax;//the largest theta that is contained within the data

        twoVector messg;//measured g(theta,phi); theta; phi
        twoVector calc; //calculated g(theta, phi);theta, phi
        twoVector grid; //grid that is used for the expansion

        /**
        *calculation of a fermi function like apodizing function
        *
        *@param gTetaPhi    intencity value g(theta, phi)
        *@param theta       value theta at gTetaPhi/g(theta, phi)
        */
        inline dataType apofct(dataType gThetaPhi, dataType theta);





    private:
        const char* infile;     //name of file with input data
        const bool apo;         //aposization falg
        const int INORM;
        std::ifstream dataF;    //for input
        std::ofstream dataW;    //for output

        /**
        * calculates the chi function from the measured data
        * according to the procedure of chuck fadley, i.e. by
        * normalizing the data on each theta-circle to zero average.
        * the Function is alway called, however it only modifies the
        * measured data if the local flag inorm < 0
        *
        */
        void calcchi();



};

#endif // PARAM_H
