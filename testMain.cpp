//#include <boost/lambda/lambda.hpp>
#include <iostream>
//#include <iterator>
//#include <algorithm>
#include <boost/math/special_functions/legendre.hpp>
#include <vector>
#include <complex>


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

    /*for(int i=0;i<10;i++)
    {

        for(int j=0;j<=i;j++)
        {
            std::cout<<"i: "<<i<<" j: "<<j;
        }
        std::cout<<std::endl;
    }*/

 /*   std::vector<double> sizing;
    std::cout<<"max size of initialized vector"<<sizing.max_size()<<std::endl;
    std::cout<<"size of initialized vector"<<sizing.size()<<std::endl;
    sizing.resize(10);
    std::cout<<"max size after resize of vector"<<sizing.max_size()<<std::endl;
    std::cout<<"size after resize of vector"<<sizing.size()<<std::endl;


    std::vector<std::vector<double> > temp;

    for(int i=0;i<10;i++)
    {
        std::vector<double> second (4,100);
        temp.push_back(second);
    }

    for(std::vector<std::vector<double> >::iterator it=temp.begin();it!=temp.end();++it)
    {
        std::cout<<(*it)[1]<<std::endl;
        /*
        for(std::vector<double>::iterator itt=(*it).begin();itt!=(*it).end();++itt)
        {
            std::cout<<*itt<<std::endl;
        }*/


    std::complex<double> cOne (2,2);
    std::cout<< cOne/2.0;
    return 0;
}
