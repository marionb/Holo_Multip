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
        *Read data file and wirite the data in to messg
        */
        void readData();

        /**
        * write the contens of messg in to output
        * @param outpur name of out put file
        */
        //void writeData(std::string output);
        static int getLineNum(std::string fileGiven);



    private:

        int const maxangles; //number of data points

        const char* infile; //name of file with input data
        std::ifstream dataF;

        double **messg; //"D dynamic array
};

#endif // PARAM_H
