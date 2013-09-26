//#include <boost/lambda/lambda.hpp>
#include <iostream>
//#include <iterator>
//#include <algorithm>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>




int main()
{


    /*for (int i=0;i<10;i++)
    {
        double x=double(i)/10.0;
        double leg=boost::math::legendre_p <double> (1,1,x);
        std::cout<< "this is the result "<<x<<"->"<<leg<<"\n";
    }

    std::vector<std::vector<double> > alm2;
    std::vector<std::vector<double> > alm1(5, std::vector<double>());
    alm2=alm1;
    for(int i=0;i<5;i++)
    {
        alm2[i].push_back(0.1);
        alm2[i].push_back(0.2);

    }
    for(int i=0;i<5;i++)
    {
        for(int j)
    }*/

    for(int i=0;i<10;i++)
    {

        for(int j=0;j<=i;j++)
        {
            std::cout<<"i: "<<i<<" j: "<<j;
        }
        std::cout<<std::endl;
    }


    return 0;
}
