//CPP Holo_Multipol
#include <iostream>
#include <cmath>
#include <complex>

#include "Multipole.h"
#include "Data.h"

int main(int argc, char** argv)
{
    std::cout<<"runing "<<argv[0];
    std::cout<<"\n----------------------------------------------------\n";


    //parameters
    std::string paramFile= "";
    std::string dataFile="";
    std::string gridData="";//contains the grid accoring to which the multipole expansion is made
    std::string inFileALM="";
    std::string outFileALM="";
    std::string outFileMPL="";
    std::string apodization=""; //apodization Flag
    int lmax=0; //max number of parameters to expand
    int symm=0; //symmetry of the tested material(or similar)
    int inorm=0;
    bool apo=false;




    if( argc == 2)
    {
        paramFile= argv[1];
    }
    else
    {
        paramFile = "param.in";
        std::cout<<"\nno parameter file given; will use "<< paramFile <<std::endl;
    }

    //read arguments from parameter File
    //TODO replace if conditions with exeptin for nicer clean up in case of error!
    std::ifstream pFile (paramFile.c_str());
    if(pFile.is_open())
    {
        std::cout<<"---------------------------------------------------------------\n";
        std::cout<<"reading file "<< paramFile <<std::endl;
        pFile>>dataFile;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
        pFile>>gridData;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
        pFile>>outFileALM;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
        pFile>>outFileMPL;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
        pFile>>lmax;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
        pFile>>symm;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
        pFile.close();
        std::cout<<"read "<< paramFile <<" successfully"<<std::endl;
        std::cout<<"--------------------------------------------------------------- \n";
    }
    else {std::cout<<"ERROR: Unable to open file: \nPlease check that it exists in the executing folder!\n";}

    //Now functions to run the computations are called.

    Multipole *file=new Multipole(dataFile, apo, inorm, lmax, symm);

    file->readData();
    file->multpl();
    file->writeAlm(outFileALM);

    if(gridData!="")
    {
        file->readGrid(gridData);
    }

    if(inFileALM!="")
    {
        file->readAlm(inFileALM);
    }

    file->expans();

    std::cout<< "write output";
    file->writeData(outFileMPL);
    //file->doyzimage(grid,0);

    //std::cout<<"apofact is "<<apofct(file->messg[0][1],file->messg[0][1]);

    delete file;
    return 0;
}








