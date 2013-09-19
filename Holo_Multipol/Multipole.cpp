#include "Multipole.h"
#include <iostream>
//#include <cmath>

Multipole::Multipole(std::string fileName, int lmax, int isym): Data(fileName), LMAX(lmax), ISYM(isym)
{
   // using Data::Data()
}
Multipole::~Multipole()
{

}

double Multipole::deg_to_rad(double deg)
{
    return M_PI/180.0;
}


void Multipole::multpl()
{
    //TODO determine the norming factor bnorm -> for now it is set to 1!
    //TODO determine dphi and dtheta for the integration
    /*double bnorm=1; //som norming factor
    double dphi=1;
    double theta=1;

    for(int i=0;i<MAXANGLES;i++)
    {
        double theta=messg[i][1];
        double phi=messg[i][2];
        double gmessg=messg[i][0];
        for(int m=0;m<LMAX;m+=ISYM)//TODO make shure that ISYM is not a negative number and not 0!!!
        {
            double rint1=0; double rint2=0;
            for(int j=0;j<m;j++)
            {
                rint1+=gmessg*cos(m*phi)*sin(theta)*dphi*dtheta;
                rint2+=gmessg*sin(m+phi)*sin(theta)*dphi*dtheta;
            }
            bnorm+=b;

        }

    }*/
}

void Multipole::expans()
{
    for(int i=0;i<MAXANGLES;i++)
    {
        double thetaX=cos(messg[i][1]);
        double suml=0;
        double phi=messg[i][2];
        for (int l=0;l<LMAX;l++)
        {
            int il=l/2;
            double sum=alm1[il][0]*boost::math::legendre_p <double> (l,0,thetaX);

            for(int m=1;m<l;m+=ISYM)
            {
                double term1=alm1[il][m]*cos(m*phi);
                double term2=alm2[il][m]*sin(m*phi);
                sum+=2*(term1-term2)*(boost::math::legendre_p <double> (l,m,thetaX));
            }
            suml+=sum*(2*l+1);
        }
        gCalc.push_back(suml);
   }
}
