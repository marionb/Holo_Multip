#include "Data.h"
#include <iostream>
//#include <iterator>

Data::Data(std::string fileName):MAXANGLES(Data::getLineNum(fileName)), infile(fileName.c_str())
{
    //the constructor takes the file name and from this name calls the static function
    //Data::getLineNum(filename) in ordet to determine the number of datapoints.
    std::cout<<"processing data file: "<<infile;
    std::cout<<" with "<<MAXANGLES<<" data points"<<std::endl;

    //initialize data arrays
    this->messg=new double[MAXANGLES][3];
    this->messgFinal=new double[MAXANGLES][3];
    thmax=0;
}

void Data::readData()
{

    dataF.open(infile);
    if(dataF.is_open())
    {
        std::cout<<"reading file";
        for(int i=0;i<MAXANGLES;i++)
        {
            //std::cout<<"read data line "<<i<<std::endl;
            dataF>>messg[i][0];
            //std::cout<<" entry 0 "<<messg[i][0];
            this->dataF>>messg[i][1];
            if(thmax<messg[i][1]) // find the maximal theta angle in the list
            {
                thmax=messg[i][1];
            }
            this->dataF>>messg[i][2];
            this->dataF.ignore (std::numeric_limits<std::streamsize>::max(), '\n'); //the rest of the row is ignored
            i++;
        }
    }else
    {
        std::cout<<"Can not open file! Will exit!";
        //Throw error here!!!
        return;

    }
    std::cout<<"read file succesfully!";
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
        std::string unused;
        while(std::getline(dataF, unused))
        {
            maxangles++;
        }
        std::cout<<fileGiven<<" has "<<maxangles<<" lines\n";

    }
    dataF.close();
    return maxangles;
}

void Data::apofct()
{
    //this function has t be caled if the certain flag is set to true!
    //what is thetamax?is it just the maximum of all thetas?
    double width =2.0; //where does this number come from?? is it always 2?

    for(int i=0;i<MAXANGLES;i++)
    {
        messg[i][0]*=1.0/(std::exp((messg[i][1]-thmax+3.0*width)/width) +1.0);
    }
}

void Data::calcchi()
{
    //TODO
    //double sumphi=0;

    //???
}

void Data::writeData(std::string output)
{
    std::cout<<"writing data into "<<output<<std::endl;

    const char* file = output.c_str();
    dataW.open(file);
    if(dataW.is_open())
    {
        dataW<<"#Intensity";
        dataW<<"    ";
        dataW<<"Theta angle";
        dataW<<"    ";
        dataW<<"Phi angle";
        dataW<<std::endl;
        for(int i=0;i<MAXANGLES;i++)
        {
            dataW<<messgFinal[i][0];
            dataW<<"    ";
            dataW<<messgFinal[i][1];
            dataW<<"    ";
            dataW<<messgFinal[i][2];
            dataW<<std::endl;
        }

    }


}
Data::~Data()
{
    delete [] messg;
    delete [] messgFinal;

}
