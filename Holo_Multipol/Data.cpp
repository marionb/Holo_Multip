#include "Data.h"
#include <iostream>
//#include <iterator>

Data::Data(std::string fileName):MAXANGLES(Data::getLineNum(fileName)), ROWNUMBER(4) ,infile(fileName.c_str())
{
    //the constructor takes the file name and from this name calls the static function
    //Data::getLineNum(filename) in ordet to determine the number of datapoints.
    std::cout<<"processing data file: "<<infile;
    std::cout<<" with "<<MAXANGLES<<" data points"<<std::endl;
    //initialize data arrays
    messg.resize( MAXANGLES, std::vector<double>( ROWNUMBER , 0 ) );
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
            //std::cout<<"\n read data line "<<i;
            dataF>>messg[i][0];
            //std::cout<<" entry 0 "<<messg[i][0];
            this->dataF>>messg[i][1];
            messg[i][1]=deg_to_rad(messg[i][1]);
            //std::cout<<" entry 1 "<<messg[i][1];
            if(thmax<messg[i][1]) // find the maximal theta angle in the list
            {
                thmax=messg[i][1];
            }
            this->dataF>>messg[i][2];
            messg[i][2]=deg_to_rad(messg[i][2]);
            //std::cout<<" entry 2 "<<messg[i][2]<<std::endl;
            this->dataF>>messg[i][3];
            messg[i][3]=deg_to_rad(messg[i][3]);
            this->dataF.ignore (std::numeric_limits<std::streamsize>::max(), '\n'); //the rest of the row is ignored
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
    //DODO Thest this function
    double sumphi=0;
    double theta=messg[0][1];
    int count=0;
    std::vector<double> normFactor;
    for (int i=0;i<MAXANGLES;i++)
    {
        sumphi+=messg[i][0];
        if(theta!=messg[i+1][1])
        {

            if(i-count<=0)
            {
              sumphi=sumphi/(i-count);
            }else{ std::cout<< "ERROR: can not calculate average in function calcchi, result migt be wrong!";}
            for(int k=count;k<=i;k++)
            {
                messg[k][0]=messg[k][0]/sumphi;
            }
            sumphi=0;
            theta=messg[i+1][1];
            count=i;
        }

    }
}

void Data::writeData(std::string output)
{
    std::cout<<"writing data into "<<output<<std::endl;

    const char* file = output.c_str();
    dataW.open(file);
    if(dataW.is_open())
    {
        dataW<<"Titel TODO\n";
        dataW<<"calc. with holo1\n";
        dataW<<"#Intensity";
        dataW<<"    ";
        dataW<<"Theta angle";
        dataW<<"    ";
        dataW<<"Phi angle\n";

        for(twoVector::iterator it=calc.begin();it!=calc.end();++it)
        {
            //std::cout<<*it[0]<<" "<<*it[1]<<" "<<*it[2]<<std::endl;
            dataW<<(*it)[0];
            dataW<<" ";
            dataW<<rad_to_deg((*it)[1]);
            dataW<<" ";
            dataW<<rad_to_deg((*it)[2]);
            dataW<<"\n";
        }

    }
}
Data::~Data()
{


}

//------------------------------------------------------//
//Static functions
double Data::deg_to_rad(double deg)
{
    return deg*M_PI/180.0;
}
double Data::rad_to_deg(double rad)
{
    return rad*180/M_PI;
}
