#include "Data.h"
#include <iostream>

Data::Data(std::string fileName):maxangles(Data::getLineNum(fileName)), infile(fileName.c_str())
{
    //the constructor takes the file name and from this name calls the static function
    //Data::getLineNum(filename) in ordet to determine the number of datapoints.
    std::cout<<"processing data file: "<<infile<<std::endl;
    std::cout<<"with "<<maxangles<<" data points"<<std::endl;
    //initialize data array
    this->messg=new double*[maxangles];
    for(int i=0;i>maxangles;i++)
    {
        this->messg[i]=new double[3];
    }

}

void Data::readData()
{

    dataF.open(infile);
    if(dataF.is_open())
    {
        std::cout<<"reading from file ...";
        for(int i=0;i<maxangles;i++)
        {
            std::cout<<"read data line "<<i<<std::endl;
            std::cout<<"data: ";
            this->dataF>>messg[i][0];
            std::cout<<messg[i][0];
            this->dataF>>messg[i][1];
            std::cout<<messg[i][1];
            this->dataF>>messg[i][2];
            std::cout<<messg[i][2];
            this->dataF.ignore (std::numeric_limits<std::streamsize>::max(), '\n'); //the rest of the row is ignored
            std::cout<<"read data line "<<i<<std::endl;
        }
    }else
    {
        std::cout<<"Can not open file input.dat! Will exit!";
        //Throw error here!!!
        return;

    }
    std::cout<<"read file succesfully! Will close file";
    this->dataF.close();
}
int Data::getLineNum(std::string fileGiven)
{
    /*
    TODO: needs exeption if file can not be opened
    */
    const char* file = fileGiven.c_str();

    int maxangles=0;
    std::ifstream dataF;
    dataF.open(file);
    if(dataF.is_open())
    {
        std::cout<<"fille open will process";
        std::string unused;
        while(std::getline(dataF, unused))
        {
            maxangles++;
        }
        std::cout<<" maxangles: "<<maxangles<<"\n";

    }
    dataF.close();
    return maxangles;
}
/*
void Data::writeData(std::string output)
{

}*/
Data::~Data()
{
    for(int i=0;i<maxangles;i++)
    {
        delete [] messg;
    }
}
