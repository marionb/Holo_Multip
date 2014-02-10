#include "Multipole.h"
#include <iostream>
#include <stdlib.h>     /* abs */

Multipole::Multipole(std::string fileName, bool apo, int inorm, int lmax, int isym): Data(fileName, apo, inorm),LMAX(lmax), ISYM(isym), INORM(inorm)//, k(2*M_PI*sqrt(ekin/150))
{
    assert(LMAX<Multipole::MAX_COEFF);
}
Multipole::~Multipole()
{
}

const int Multipole::getLMAX()
{
    return LMAX;
}

void Multipole::multpl()
{
    std::cout<<"\n--------------------------------------------------------------- \ncalculating multipole coefficients"<<std::endl;
    dataType bnorm = 1;
    alm.clear();
    for(int l=0;l<=LMAX;l+=2)
    {
        std::vector<std::complex<dataType> > lcoeff;
        for(int m=0;m<=l;m+=ISYM)
        {
            std::complex<dataType> int_res=0;


            //go through all angles alm=int(g(theta,phi)*conj(Y_lm(thea,phi)) dOmega
            for(int i=0;i<MAXANGLES;i++)
            {
                dataType phi=messg[i][2];
                dataType theta=messg[i][1];
                dataType gmessg=messg[i][0];
                dataType dOmega=messg[i][3];

                std::complex<float> temp1=std::conj(boost::math::spherical_harmonic<dataType>(l, m, theta, phi))*dOmega;
                std::complex<float> temp=gmessg*temp1;
                int_res+=temp;

            }

            if(l==0 && m==0 && std::real(int_res)>0.001)
            {
                bnorm=std::real(int_res);
                std::cout<<"bnorm="<<bnorm<<std::endl;
            }
            lcoeff.push_back(int_res/bnorm);

        }

        alm.push_back(lcoeff);
    }
    removeAlmBackground();
}

void Multipole::expans()
{
    std::cout<<"\n--------------------------------------------------------------- \nexpanding function according to calculated coefficients \n";
    calc.clear();

    std::vector< std::vector<dataType> >::const_iterator itr;
    for (itr = grid.begin(); itr != grid.end(); ++itr)
    {
        dataType theta=(*itr)[0];
        dataType phi=(*itr)[1];
        dataType g_theta_phi=intencity(theta,phi);

        std::vector<dataType> temp;

        temp.push_back(g_theta_phi);
        temp.push_back(theta);
        temp.push_back(phi);

        calc.push_back(temp);
   }
}

dataType Multipole::intencity(dataType theta, dataType phi)
{
    std::complex<dataType> summ;//=(0,0);
    dataType summ0 = 0;
    for(int l=0; l<=LMAX; l+=2)
    {
        int m=0;
        std::complex<dataType> coeff=alm[l/2][m];
        std::complex<dataType> ylm=boost::math::spherical_harmonic<dataType> (l,m,theta,phi);
        summ0+=std::real(coeff*ylm);
        m+=ISYM;
        for(;m<=l;m+=ISYM)
        {
            summ+=alm[l/2][m]*boost::math::spherical_harmonic<dataType> (l,m,theta,phi);
        }
    }

    dataType totSumm=summ0+2*std::real(summ);
    return totSumm;
}

void Multipole::removeAlmBackground()
{
    if(INORM<0)
    {
        //Function does nothing if the flag INORM is less than 0
        return;
    }
    std::cout<<"removing alm with l<"<<INORM<<std::endl;
    for(int i=0;i<=INORM;i+=2)
    {
        int il=i/2;
        for(int j=0;j<=i;j++)
        {

            alm[il][j]=0;

        }
    }

}

/*READ WRITE functions*/
void Multipole::writeAlm(std::string almFile)
{

    std::cout<<"\n--------------------------------------------------------------- \n";
    std::cout<<"writing A_lm to file "<<almFile<<std::endl;

    const char* file = almFile.c_str();
    std::ofstream almWrite;
    almWrite.open(file);
    if(almWrite.is_open())
    {
        almWrite<<"#Coeficients for multipole expansion A_lm\n";
        almWrite<<"#Calculated using Holo_Multipol\n";
        almWrite<<"#Using ";
        almWrite<<LMAX;
        almWrite<<" coefficients\n";
        almWrite<<"#l";
        almWrite<<"    ";
        almWrite<<"m";
        almWrite<<"    ";
        almWrite<<"A_lm\n";
        for(int i=0;i<=LMAX;i+=2)
        {
            int il=i/2;
            for(int j=0;j<=i;j++)
            {
                almWrite<<i;
                almWrite<<" ";
                almWrite<<j;
                almWrite<<" ";
                almWrite<<alm[il][j];
                almWrite<<"\n";
            }
        }

    }
    else
    {
        std::cout<<"can not write coefficients to file! "<<almFile<<std::endl;
        return;
    }
    std::cout<<"coefficients for multipole expansion can be found in "<<almFile<<std::endl;
    std::cout<<"--------------------------------------------------------------- \n";
}

void Multipole::readAlm(std::string almFile)
{

    std::cout<<"\n--------------------------------------------------------------- \n";
    std::cout<<"readig A_lm from "<<almFile<<std::endl;

    const char* file = almFile.c_str();
    std::ifstream almRead;
    almRead.open(file);
    if(almRead.is_open())
    {
        alm.clear();
        std::string line;
        std::vector<std::complex<dataType> > lcoeff;
        char dummy='#';
        int lOld,l,m;
        lOld=0;
        l=0;
        m=0;
        std::complex<dataType> coeff;
        while(std::getline(almRead,line))
        {
            //ignore all lines that contain a #
            if(line.find(dummy)!=std::string::npos)
            {
                continue;
            }

            std::istringstream ss(line);

            ss>>l>>m>>coeff;

            //if l is increast append the temporary vector to the alm vector
            if(lOld!=l)
            {
                alm.push_back(lcoeff);
                lcoeff.clear();
            }

            lcoeff.push_back(coeff);

            lOld=l;

        }

    }
    else
    {
        std::cout<<"can not read coefficients from file! "<<almFile<<std::endl;
        return;
    }
    std::cout<<"--------------------------------------------------------------- \n";
}
