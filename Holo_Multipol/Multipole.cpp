#include "Multipole.h"
#include <iostream>
#include <stdlib.h>     /* abs */

Multipole::Multipole(std::string fileName, int lmax, int isym, double ekin): Data(fileName),LMAX(lmax), ISYM(isym), k(2*M_PI*sqrt(ekin/150))
{
    assert(LMAX<Multipole::MAX_COEFF);
    alm1.resize( LMAX/2+1, std::vector<double>());
    alm2.resize( LMAX/2+1, std::vector<double>());
}
Multipole::~Multipole()
{}

void Multipole::multpl()
{
    //TODO: spherical_harmonic needs to be calculated as spherical_harmonic* (complex conjugated)
    double omega=4.0*M_PI;

    std::cout<<"\n in multipl function"<<std::endl;

    for(int l=0;l<=LMAX;l+=2)
    {
        int il=l/2;

//        std::cout<<"(l, il)"<<l<<","<<il<<std::endl;

        //alm1[il].push_back(std::vector<double>(il)); //append a vector of length l/2+1 to make shure there is at least one element; initialize with zeros
        //alm2[il].push_back(std::vector<double>(il));
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
            alm2[il].push_back(rint2/MAXANGLES*omega);
            //std::cout<<"result final"<<alm1[il][m];
        }
        //std::cout<<std::endl;
    }


}

void Multipole::expans()
{
    std::cout<<"in expans function\n";
    for(int i=0;i<MAXANGLES;i++)
    {
        std::cout<<"i "<<i<<std::endl;
        double theta=messg[i][1];
        double phi=messg[i][2];
        double suml=0;
        std::cout<<"(theta, phi)=("<<theta<<", "<<phi<<")\n";
        for (int l=0;l<LMAX;l+=2)
        {
            //std::cout<<"l is "<<l<<std::endl;
            int il=l/2;
            double sum=alm1[il][0]*boost::math::spherical_harmonic_r <double> (l,0,theta,phi); //m=0
            //std::cout<<"sum"<<sum<<std::endl;
            for(int m=1;m<=l;m++)
            {
                std::cout<<"in seccond for loop!";
                //calculate RE(sum_m A_lm*Y_lm)=sum_m RE(A_lm)*RE(Y_lm)-IM(A_lm)*IM(Y_lm)
                double real=alm1[il][m]*boost::math::spherical_harmonic_r <double> (l,m,theta, phi); //imaginary part
                double imag=alm2[il][m]*boost::math::spherical_harmonic_i <double> (l,m,theta, phi); //real part
                sum+=2*(real-imag); //
                //std::cout<<"sum, i"<<sum<<", "<<l<<std::endl;
            }
        }
        calc[i][0]=suml;
        calc[i][1]=suml;
        calc[i][2]=suml;
   }
}


void Multipole::holorad(double alpha, double beta)
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
    }*/
}


/*-------------------Private--------------------------------------*/
double Multipole::innerSumm(int l, double alpha, double beta)
{
    double summ=0;
    int il=l/2;
    for(int m=ISYM;m<=l;m+=ISYM)
        {
            double real=alm1[il][m]*boost::math::spherical_harmonic_r <double> (l,m,alpha, beta); //imaginary part
            double imag=alm2[il][m]*boost::math::spherical_harmonic_i <double> (l,m,alpha, beta);
            summ+=2*(real-imag);
        }
    return summ;
}

void Multipole::printAlm()
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
}
