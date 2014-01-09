#include "Multipole.h"
#include <iostream>
#include <stdlib.h>     /* abs */

Multipole::Multipole(std::string fileName, int lmax, int isym, double ekin): Data(fileName),LMAX(lmax), ISYM(isym), k(2*M_PI*sqrt(ekin/150))
{
    assert(LMAX<Multipole::MAX_COEFF);
    alm.resize( LMAX/2+1, std::vector<std::complex<double> >());
}
Multipole::~Multipole()
{
    //TODO enentiuelle leichen löschen
    //delete alm1;
    //delete alm2;
}

const int Multipole::getLMAX()
{
    return LMAX;
}

void Multipole::multpl()
{
    std::cout<<"\n --------------------------------------------------------------- \n expanding multipole coefficients"<<std::endl;
    double bnorm = 1;
    for(int l=0;l<=LMAX;l+=2)
    {
        int il=l/2;

//        std::cout<<"(l, il)"<<l<<","<<il<<std::endl;

        //alm1[il].push_back(std::vector<double>(il)); //append a vector of length l/2+1 to make shure there is at least one element; initialize with zeros
        //alm2[il].push_back(std::vector<double>(il));
        for(int m=0;m<=l;m+=ISYM)
        {
            std::complex<double> int_res=0;

            //go through all angles alm=int(g(theta,phi)*conj(Y_lm(thea,phi)) dOmega
            for(int i=0;i<MAXANGLES;i++)
            {
                double phi=messg[i][2];
                double theta=messg[i][1];
                double gmessg=messg[i][0];
                double dOmega=messg[i][3];

                //std::cout<<"i, phi, theta, g "<<i<<" "<<phi<<","<<theta<<", "<<gmessg<<std::endl;
                int_res+=gmessg*std::conj(boost::math::spherical_harmonic<double>(l, m, theta, phi))*dOmega;


                //std::cout<<"real, imag "<<rint1<<", "<<rint2<<"\n";


            }
            if(l==0 && m==0 && std::real(int_res)>0.001)
            {
                bnorm=std::real(int_res);
                std::cout<<"bnorm="<<bnorm<<std::endl;
            }
            alm[il].push_back(int_res/ bnorm);
            std::cout<<l<<" "<<m<<" "<<alm[il][m]<<std::endl;//<<" "<<alm2[il][m]<<std::endl;
        }
        //std::cout<<std::endl;
    }
std::cout<<"alm: (l,m,real,imag)\n";

}

void Multipole::expans()
{
    std::cout<<"\n --------------------------------------------------------------- \n expanding function according to calculated coefficients \n";
    calc.clear();
    for(int i=0;i<MAXANGLES;i++)
    {
        std::cout<<"i "<<i<<std::endl;
        double theta=messg[i][1];
        double phi=messg[i][2];
        double g_theta_phi=intencity(theta,phi);
        std::vector<double> temp;
        std::cout<<g_theta_phi<<" "<<180/M_PI*theta<<" "<<180/M_PI*phi<<std::endl;
        temp.push_back(g_theta_phi);
        temp.push_back(theta);
        temp.push_back(phi);//TODO if time: make this call better -> works for a start

        calc.push_back(temp);
   }
   std::cout<<"calc:(g(theat,phi),theta,phi)\n";
}


/*void Multipole::holorad(double alpha, double beta)
{
    double r=0;
    double kr=0;
    double suml=0;
    for(int nr=1;nr<=100;nr++)
    {
        r=nr*0.1;
        kr=k*r;
        for(int l=0;l<LMAX;l+=2)
        {
            int il=l/2;
            suml+=innerSumm(l,alpha,beta)*boost::math::sph_bessel (l, kr)*vorz(il); //TODO:not shure if the order of l is correct
        }
        radimg.push_back(suml);
    }

}

void Multipole::doyzimage(double grid, int xyz)
{
    double x=xyz*grid; //=0!
    double alpha=0;
    double beta=0;
    double suml=0;
    int max=63;

    image2D.resize(max);
    for(int nz=0;nz<=max;nz++)
    {
        image2D[nz].expand(max);
        double z=nz*grid;
        for(int ny=-max;ny<=max;ny++)
        {
            double y=ny*grid;
            double rho=sqrt(x*x+y*y);
            double rVec=sqrt(rho*rho+z*z);
            alpha=calcth(rho, z);
            beta=calcphi(x,y);
            double kr=k*rVec;
            for(int l=0;l<=LMAX;l+=2)
            {
                int il=l/2;
                suml+=innerSumm(l,alpha,beta)*boost::math::sph_bessel (l, kr)*vorz(il);
            }
            (ny<0)?(image2D[nz].negative[abs(ny)]=suml):(image2D[nz].negative[ny]=suml); //if ny<0 safe in negative part else (if ny>=0 safe in positive part)
        }
    }
}

void Multipole::doxyzimage(double grid)
{
    //TODO test if everything works fine with vector structure
    int max=63;
    image3D.resize(max);
    for(int nz=0;nz<=max;nz++)
    {
        double z=nz*grid;
        image3D[nz].expand(max);
        for(int ny=-max;ny<=max;ny++)
        {
            double y=ny*grid;
            //(ny<0)? (image3D[nz].negative[-ny].expand(max)):(image3D[nz].positive[ny].expand(max)); //expand only the enty ny that is used next
            SpecialVector<double> xtemp;
            xtemp.expand(max);
            for(int nx=-max;nx<=max;nx++)
            {
                double x=nx*grid;
                double r=sqrt(x*x+y*y+z*z);
                double alpha=acos(z/r);
                double beta=calcphi(x,y);
                double kr =k*r;
                double suml=0;
                for(int l=0;l<=LMAX;l+=2)
                {
                    int il=l/2;
                    suml+=innerSumm(l,alpha,beta)*boost::math::sph_bessel (l, kr)*vorz(il);
                }
                xtemp.addElement(suml,nx);
            }
            image3D[nz].addElement(xtemp,ny);
        }
    }
    image3D[0].positive[0].addElement(0,0);
}


void Multipole::scaleimage(double grid)
{
    int max=63;
    assert(image3D.size()!=0);
    for(int iz=0;iz<=max;iz++)
    {
        for(int iy=-max;iy<max;iy++)
        {
            SpecialVector<double> xtemp;
            xtemp.expand(max);
            for(int ix=-max;ix<max;ix++)
            {
                double r=grid*sqrt(double(ix*ix+iy*iy+iz*iz));
                double temp= r*image3D[iz].getElement(iy).getElement(ix);
                xtemp.addElement(temp, ix);

            }
            image3D[iz].addElement(xtemp,iy);
        }
    }

}

void Multipole::smooth(double grid)
{
   /* assert(image3D.size()!=0);
    int box = 7; //smoothing box size half
    double width = 1; //gauss sigma width
    std::cout<<"pewforming gaussian smoothing\n";
    std::vector<SpecialVector<SpecialVector<double> > > >simage;
    simage=image3D;
    for(int nz=0;nz<=63;nz++)
    {
        for(int ny=0;ny<=63;ny++)
        {
            for(int nx=0;nx<=63;nx++)
            {
                double imageR=image3D[nz].getElement(ny).getElement(nx);
                for(iz=-box;iz<=box;iz++)
                {
                    for(iy=-box;iy<=box;iy++)
                    {
                        for(ix=-box;ix<=box;ix++)
                        {
                            double distsq=(ix+iy+iz)*(ix+iy+iz);
                            double gs=1/(width*2.5)*exp(-0.5*distsq/(width*width));
                            imageR+=sigma
                        }
                    }
                }
            }
        }
    }*--->/
}


/*-------------------Private--------------------------------------*/
double Multipole::intencity(double theta, double phi)
{
    double summ0=0;
    std::complex<double> summ (0,0);
    for(int l=0; l<=LMAX; l+=2)
    {
        int m=0;
        summ0+=std::real(alm[l/2][m]*boost::math::spherical_harmonic<double> (l,m,theta,phi));
        m+=ISYM;
        for(;m<=l;m+=ISYM)
        {
            summ+=alm[l/2][m]*boost::math::spherical_harmonic<double> (l,m,theta,phi);
        }
    }


    return summ0+2*std::real(summ);
}

/*void Multipole::printAlm()
{
    for(unsigned int i=0;i<alm1.max_size();i++)
    {
        for(unsigned int j=0;j<alm1[i].max_size();j++)
        {
            std::cout<<alm1[i][j]<<";";
        }
        std::cout<<"i \n";
    }
}

int Multipole::vorz(int exp)
{
    exp=abs(exp);
    if(exp%2==0)
    {
        return 1;
    }else return -1;
}

double Multipole::calcth(double y, double z)
{
    if(y==0)
    {
        return 0;
    }else
    {
        if(z<0)
        {
            return M_PI-abs(atan(abs(y/z)));
        }else if(z==0.0)
        {
            return M_PI;
        }else if(z>0)
        {
            return abs(atan(abs(y/z)));
        }
        return -1;
    }
}

double Multipole::calcphi(double x, double y)
{
    if(x<0)
    {
        if(y<0)
        {
            return M_PI+abs(atan(y/x));
        }else if(y==0.0)
        {
            return M_PI;
        }else if(y>0)
        {
            return M_PI-abs(atan(y/x));
        }
    }else if(x==0.0)
    {
        if(y<0)
        {
            return 3.0/2.0*M_PI;
        }else if(y==0.0)
        {
            return 0;
        }else if(y>0)
        {
            return M_PI;
        }
    }else if(x>0)
    {
        if(y<0)
        {
            return 2.0*M_PI-abs(atan(y/x));
        }else if(y==0.0)
        {
            return 0;
        }else if(y>0)
        {
           return atan(y/x);
        }
    }
    return -1;
}*/
