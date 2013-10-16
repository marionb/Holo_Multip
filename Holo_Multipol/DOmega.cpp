#include "DOmega.h"

DOmega::DOmega(std::string file):gridFile(file.c_str())
{

}

void DOmega::readFile()
{
    grid.open(gridFile);
    if(grid.is_open())
    {
        std::string line;
        while (std::getline(grid, line))
        {
            std::string str ("#");
            if(line.find(str)!=std::string::npos) continue; //all lines containing a "#" are thought to be coments and are not read

            std::istringstream iss(line);
            int index;
            double theta, dePhi;

            if (!(iss >> index >> theta >> dePhi))
            {
                std::cout<<"error reading";
                break;
            } // error

            Triple <int, double, double> temp;
            temp.first=index;
            temp.second=theta;
            temp.third=dePhi;
            dphi.push_back(temp);
        }
    }else
    {
        std::cout<<"Can not open file! Will exit!";
        //Throw error here!!!
        return;

    }
    std::cout<<"read file succesfully!";
    grid.close();
}


void DOmega::set_dtheta()
{
    if(dphi.empty())
    {
        std::cout<<"ERROR: no values available to calculate dtheta!";
        return;
    }
    dtheta=abs(dphi[0].second-dphi[1].second);
}

void DOmega::fitToGrid(Data element)
{
    for(int i=0;i<element.MAXANGLES;i++)
    {

        element.messg[i][]
    }

}
