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

    private:

        const int MAXANGLES; //number of data points

        const char* infile; //name of file with input data
        std::ifstream dataF; //for input
        std::ofstream dataW; //for output

    protected:

        double (*messg)[3]; //dynamic array containing input data -> this data needs to be accesible
        double (*messgFinal)[3]; //dynamic array containing output data
        double thmax;
};

#endif // PARAM_H
