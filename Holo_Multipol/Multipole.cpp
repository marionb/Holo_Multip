#include "Multipole.h"
#include <iostream>
#include <stdlib.h>     /* abs */

Multipole::Multipole(std::string fileName, int lmax, int isym, double ekin): Data(fileName),LMAX(lmax), ISYM(isym)
{
    k=2*M_PI*sqrt(ekin/150);
    alm1.resize( MAXANGLES, std::vector<double>());
    alm2.resize( MAXANGLES, std::vector<double>());
}
Multipole::~Multipole()
{

}

std::vector<double> Multipole::get_radimg()
{
    return radimg;
}
void Multipole::multpl()
{
    //TODO: spherical_harmonic needs to be calculated as spherical_harmonic* (complex conjugated)
    double omega=4.0*M_PI;

    std::cout<<"\n in multipl function"<<std::endl;

    for(int l=0;l<LMAX;l+=2)
    {
        int il=l/2;

        std::cout<<"(l, il)"<<l<<","<<il<<std::endl;

        alm1.push_back(std::vector<double>( il , 0 ) ); //append a vector of length l/2+1 to make shure there is at least one element; initialize with zeros
        alm2.push_back(std::vector<double>( il , 0 ) );
        for(int m=0;m<=l;m++)
        {
            double rint1=0;
            double rint2=0;
            for(int i=0;i<MAXANGLES;i++)
            {
                //To optimize do a vector of functions containing the spherical harmonics for a given l and m.
                double phi=messg[i][2];
                double theta=messg[i][1];
                double gmessg=messg[i][0];

                //std::cout<<"i, phi, theta, g "<<i<<" "<<phi<<","<<theta<<", "<<gmessg<<std::endl;
                rint1+=vorz(m)*gmessg*boost::math::spherical_harmonic_r<double>(l, (-1)*m, theta, phi);

                rint2+=vorz(m)*gmessg*boost::math::spherical_harmonic_i<double>(l, (-1)*m, theta, phi);
                //std::cout<<"real, imag "<<rint1<<", "<<rint2<<"\n";

            }
            //std::cout<<"ended inner loop will write result \n";

            alm1[il].push_back(rint1/MAXANGLES*omega);
            //alm2[il][m]=rint2/omega;
            std::cout<<"result"<<alm1[il][m]<<std::endl;
        }
    }


}

void Multipole::expans()
{
    for(int i=0;i<MAXANGLES;i++)
    {
        double theta=messg[i][1];
        double phi=messg[i][2];
        double suml=0;
        for (int l=0;l<LMAX;l+=2)
        {
            int il=l/2;
            double sum=alm1[il][0]*boost::math::spherical_harmonic_r <double> (l,0,theta,phi); //m=0

            for(int m=1;m<=l;m++)
            {
                //calculate RE(sum_m A_lm*Y_lm)=sum_m RE(A_lm)*RE(Y_lm)-IM(A_lm)*IM(Y_lm)
                double real=alm1[il][m]*boost::math::spherical_harmonic_r <double> (l,m,theta, phi); //imaginary part
                double imag=alm2[il][m]*boost::math::spherical_harmonic_i <double> (l,m,theta, phi); //real part
                sum+=2*(real-imag); //
            }
        }
        calc[i][0]=suml;
        calc[i][1]=suml;
        calc[i][2]=suml;
   }
}


void Multipole::holorad(double alpha, double beta)
{
    /*std::vector<double> summ;
    double var=cos(alpha);
    for(int l=0;l<LMAX;l+=2)
    {
        int il=l/2;
        summ.push_back(alm1[0][il]);plegndr(l,0,var));
        for(int m=ISYM;m<l;m+=ISYM)
        {
            double real=alm1[il][m]*boost::math::spherical_harmonic_r <double> (l,m,alpha, beta); //imaginary part
            double imag=alm2[il][m]*boost::math::spherical_harmonic_i <double> (l,m,alpha, beta);
            summ[il]+=2*(real-imag);
        }

    }
    double r=0;
    double kr=0;
    double suml=0;
    for(int nr=1;nr<100;nr++)
    {
        r=nr*0.1;
        kr=k*r;
        kr=k*r;
        suml=0;
        for(int l=0;l<LMAX;l+=2)
        {
            int il=l/2;
            suml+=summ(il)*boost::math::sph_bessel <double> (l, kr)*vorz(il); //TODO:not shure if this order is correct
        }
        radimg.push_back(suml);
    }*/

}

void Multipole::helperFunc(double alpha, double beta)
{

}

int Multipole::vorz(int exp)
{
    exp=abs(exp);
    if(exp%2==0)
    {
        return 1;
    }else return -1;
}
