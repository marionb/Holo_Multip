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
    std::string outFileALM="";
    std::string outFileMPL="";
    int lmax=0; //max number of parameters to expand
    int symm=0; //symmetry of the tested material(or similar)



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


    Multipole *file=new Multipole(dataFile,lmax,symm);
    file->readData();
    file->multpl();
    file->writeAlm(outFileALM);
    file->expans();

    std::cout<< "write output";
    file->writeData(outFileMPL);
    //file->doyzimage(grid,0);

    //std::cout<<"apofact is "<<apofct(file->messg[0][1],file->messg[0][1]);

    delete file;
    return 0;
}








