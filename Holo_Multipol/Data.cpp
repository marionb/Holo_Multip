#include "Data.h"
#include <iostream>
//#include <iterator>

Data::Data(std::string fileName, bool apo, bool inorm):MAXANGLES(Data::getLineNum(fileName)), ROWNUMBER(4) ,infile(fileName.c_str()),apo(apo), INORM(-inorm)
{
    //the constructor takes the file name and from this name calls the static function
    //Data::getLineNum(filename) in ordet to determine the number of datapoints.
    std::cout<<"processing data file: "<<infile;
    std::cout<<" with "<<MAXANGLES<<" data points"<<std::endl;
    //initialize data arrays
    messg.resize( MAXANGLES, std::vector<dataType>( ROWNUMBER , 0 ) );
    thmax=0;
}

void Data::readData()
{

    dataF.open(infile);
    if(dataF.is_open())
    {
        grid.clear();
        std::cout<<"--------------------------------------------------------------- \nreading file "<<infile<<std::endl;
        for(int i=0;i<MAXANGLES;i++)
        {
            dataType inNum0; // g(theta,phi)
            dataType inNum1; // theta
            dataType inNum2; //phi
            dataType inNum3; //DOmega
            //std::cout<<"\n read data line "<<i<<"   ";
            this->dataF>>inNum0;
            //std::cout<<inNum0<<" ";
            this->dataF>>inNum1;
            //std::cout<<inNum1<<" ";
            this->dataF>>inNum2;
            //std::cout<<inNum2<<" ";
            this->dataF>>inNum3;
            //std::cout<<inNum3<<" ";

            messg[i][0]=inNum0;

            //std::cout<<" entry 0 "<<messg[i][0];

            //std::cout<<" || entry theta "<<inNum1;
            //std::cout<<" entry phi "<<inNum2;
            if(thmax<inNum1) // find the maximal theta angle in the list
            {
                thmax=inNum1;
            }
            //convert angles from degres to radian

            messg[i][1]=deg_to_rad(inNum1);

            messg[i][2]=deg_to_rad(inNum2);

            messg[i][3]=deg_to_rad(inNum3);
            thmax=deg_to_rad(thmax);

            if(apo)
            {
                messg[i][0]= apofct(messg[i][0], messg[i][1]);
            }

            //fill grid vector
            std::vector<dataType> temp;
            temp.push_back(deg_to_rad(inNum1));
            temp.push_back(deg_to_rad(inNum2));
            grid.push_back(temp);



            this->dataF.ignore (std::numeric_limits<std::streamsize>::max(), '\n'); //the rest of the row is ignored
        }
    }else
    {
        std::cout<<"Can not open file! Will exit! \n---------------------------------------------------------------";
        //Throw error here!!!
        return;

    }
    std::cout<<"read file succesfully! \n---------------------------------------------------------------\n";
    this->dataF.close();
    calcchi();

}

int Data::getLineNum(std::string fileGiven)
{

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

void Data::readGrid(std::string gridFile)
{
    std::cout<<"---------------------------------------------------\n";
    std::cout<<"reading file containing grid "<<gridFile<<std::endl;

    std::ifstream grData (gridFile.c_str());
    dataType theta, phi;
    if(grData.is_open())
    {
        grid.clear(); //make shure grid is empty
        std::cout<<"grit used in expansion\n";
        while(grData>>theta >> phi)
        {
            //std::cout<< theta<<"    "<<phi<<std::endl;
            std::vector<dataType> temp;
            temp.push_back(theta);
            temp.push_back(phi);
            grid.push_back(temp);
        }
        grData.close();
    }
    else
    {
        std::cout<<"unable to open "<<gridFile<<std::endl;
        std::cout<<"using grid from input data";
        return;
    }

    std::cout<<"successfully read grid file\n";
    std::cout<<"---------------------------------------------------";
}

void Data::printGrid()
{
    std::cout<<"----------------------\n";
    std::cout<<"print grid contence\n";
    std::cout<<"theta   phi\n";
    if(grid.empty())
    {
        std::cout<<"grid is empty!\n";
        return;
    }

    std::vector< std::vector<dataType> >::const_iterator row;
    for (row = grid.begin(); row != grid.end(); ++row)
    {
        std::cout<<(*row)[0]<<"     "<<(*row)[1]<<std::endl;
    }

    std::cout<<"----------------------\n";
}

dataType Data::apofct(dataType gThetaPhi, dataType theta)
{
    //what is thetamax?is it just the maximum of all thetas?
    dataType width =2.0; //where does this number come from?? is it always 2?

    return gThetaPhi*1.0/(std::exp((theta-this->thmax+3.0*width)/width) +1.0);
}

void Data::calcchi()
{
    //TODO Test this function
    //funcion is only then called if normalization flag is set to <0
    if(INORM<0)
    {
        return;
    }
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"calculating the chi function on the measured data"<<std::endl;
    dataType sumphi=0;
    dataType theta=messg[0][1]; //set first theta value
    int count=0;


    for (int i=0;i<MAXANGLES;i++)
    {
        ++count;
        sumphi+=messg[i][0]; //add g(theta,phi)
        if((((i+1)<MAXANGLES)&&(theta!=messg[i+1][1]))or i==MAXANGLES)
        {

            if(count>0)
            {
              sumphi=sumphi/(count);
            }
            else{ std::cout<< "ERROR: can not calculate average in function calcchi, result migt be wrong!\n";}

            for(int k=count;k<=i;k++)
            {
                messg[k][0]=messg[k][0]/sumphi;
            }
            if(i<MAXANGLES)
            {
                sumphi=0;
                theta=messg[i+1][1];
                count=0;
            }

        }

    }
    std::cout<<"---------------------------------------------------------------"<<std::endl;
}

void Data::writeData(std::string output)
{
    std::cout<<"writing data into "<<output<<std::endl;

    const char* file = output.c_str();
    dataW.open(file);
    if(dataW.is_open())
    {
        dataW<<"#Expanded function\n";
        dataW<<"#Calculated using Holo_Multipol\n";
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
dataType Data::deg_to_rad(dataType deg)
{
    return deg*M_PI/180.0;
}
dataType Data::rad_to_deg(dataType rad)
{
    return rad*180/M_PI;
}
