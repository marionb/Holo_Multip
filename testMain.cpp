#include <iostream>
#include <fstream>

int getLineNum(std::string fileGiven)
{
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

int main()
{
    getLineNum("input.dat");
    return 0;
}
