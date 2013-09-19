#ifndef DATA_H
#define DATA_H

#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>

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

        friend class Multipole;


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


//        bool iapo;
        double thmax;

//    protected:

        double (*messg)[3]; //dynamic array containing input data -> this data needs to be accesible
        double (*messgFinal)[3]; //dynamic array containing output data

    private:
        const char* infile; //name of file with input data
        std::ifstream dataF; //for input
        std::ofstream dataW; //for output
};

#endif // PARAM_H
