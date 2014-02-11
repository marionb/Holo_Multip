//CPP Holo_Multipol
#include <iostream>
#include <cmath>
#include <complex>

#include "Multipole.h"
#include "Data.h"


struct ReadInput{
    //parameters used
    std::string paramFile;
    std::string dataFile;
    std::string gridData;//contains the grid accoring to which the multipole expansion is made
    std::string inFileALM;
    std::string outFileALM;
    std::string outFileMPL;
    std::string apodization; //apodization Flag
    int lmax; //max number of parameters to expand
    int symm; //symmetry of the tested material(or similar)
    int inorm;
    bool apo;

    void readfile(std::string paramFile)
    {
        std::ifstream pFile (paramFile.c_str());
        if(pFile.is_open())
        {
            std::cout<<"---------------------------------------------------------------\n";
            std::cout<<"reading file "<< paramFile <<std::endl;
            std::string dummy;

            std::getline(pFile,dummy); //file format version not needed

            pFile>>dataFile; //data input fiel (angle set)
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            std::getline(pFile,dummy); //title not needed

            pFile>>outFileALM; //output file for ALM coefficients
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            pFile>>outFileMPL; //output file for multipole
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            std::getline(pFile,dummy); //file for xy output not needed

            std::getline(pFile,dummy); //file for xyz output not needed

            std::getline(pFile,dummy); //file for radial output not needed

            pFile>>lmax; //lmax of multipole expansion
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            pFile>>symm; //symmetry
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            pFile>>inorm; //normalization flag
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            pFile>>apo; //apodization flag
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            pFile>>inFileALM;
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            std::getline(pFile,dummy);
            std::getline(pFile,dummy);
            std::getline(pFile,dummy);
            std::getline(pFile,dummy);
            std::getline(pFile,dummy);
            std::getline(pFile,dummy);
            std::getline(pFile,dummy);
            std::getline(pFile,dummy);

            pFile>>gridData;
            pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

            pFile.close();
            std::cout<<"read "<< paramFile <<" successfully"<<std::endl;
            std::cout<<"--------------------------------------------------------------- \n";
        }
        else {std::cout<<"ERROR: Unable to open file: \nPlease check that it exists in the executing folder!\n";}
    }
}input;



int main(int argc, char** argv)
{
    std::cout<<"runing "<<argv[0];
    std::cout<<"\n----------------------------------------------------\n";
    std::string dummy;
    std::string paramFile;

    if( argc == 2)
    {
        paramFile= argv[1];
    }
    else
    {
        paramFile = "input.par";
        std::cout<<"\nno parameter file given; will use "<< paramFile <<std::endl;
    }

    input.readfile(paramFile);

    //read arguments from parameter File
   /* std::ifstream pFile (paramFile.c_str());
    if(pFile.is_open())
    {
        std::cout<<"---------------------------------------------------------------\n";
        std::cout<<"reading file "<< paramFile <<std::endl;
        pFile>>dummy; //file format version not needed
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.dataFile; //data input fiel (angle set)
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>dummy; //title not needed
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.outFileALM; //output file for ALM coefficients
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.outFileMPL; //output file for multipole
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>dummy; //file for xy output not needed
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>dummy; //file for xyz output not needed
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>dummy; //file for radial output not needed
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.lmax; //lmax of multipole expansion
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.symm; //symmetry
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.inorm; //normalization flag
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.apo; //apodization flag
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>input.inFileALM;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile>>gridData;
        pFile.ignore (std::numeric_limits<std::streamsize>::max(), '\n');

        pFile.close();
        std::cout<<"read "<< paramFile <<" successfully"<<std::endl;
        std::cout<<"--------------------------------------------------------------- \n";
    }
    else {std::cout<<"ERROR: Unable to open file: \nPlease check that it exists in the executing folder!\n";}
*/
    //Now functions to run the computations are called.

    Multipole *file=new Multipole(input.dataFile, input.apo, input.inorm, input.lmax, input.symm);

    file->readData();
    file->multpl();
    file->writeAlm(input.outFileALM);

    if(input.gridData!="")
    {
        file->readGrid(input.gridData);
    }

    if(input.inFileALM!="")
    {
        file->readAlm(input.inFileALM);
    }

    file->expans();

    std::cout<< "write output";
    file->writeData(input.outFileMPL);

    delete file;
    return 0;
}








