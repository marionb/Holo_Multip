#include "Multipole.h"
#include <iostream>
#include <stdlib.h>     /* abs */

Multipole::Multipole(std::string fileName, bool apo, int inorm, int lmax, int isym): Data(fileName, apo),INORM(inorm),LMAX(lmax), ISYM(isym)//, k(2*M_PI*sqrt(ekin/150))
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
    std::cout<<"\n --------------------------------------------------------------- \n calculating multipole coefficients"<<std::endl;
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
    std::cout<<"\n --------------------------------------------------------------- \n expanding function according to calculated coefficients \n";
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


/*void Multipole::holorad(dataType alpha, dataType beta)
{
    dataType r=0;
    dataType kr=0;
    dataType suml=0;
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

void Multipole::doyzimage(dataType grid, int xyz)
{
    dataType x=xyz*grid; //=0!
    dataType alpha=0;
    dataType beta=0;
    dataType suml=0;
    int max=63;

    image2D.resize(max);
    for(int nz=0;nz<=max;nz++)
    {
        image2D[nz].expand(max);
        dataType z=nz*grid;
        for(int ny=-max;ny<=max;ny++)
        {
            dataType y=ny*grid;
            dataType rho=sqrt(x*x+y*y);
            dataType rVec=sqrt(rho*rho+z*z);
            alpha=calcth(rho, z);
            beta=calcphi(x,y);
            dataType kr=k*rVec;
            for(int l=0;l<=LMAX;l+=2)
            {
                int il=l/2;
                suml+=innerSumm(l,alpha,beta)*boost::math::sph_bessel (l, kr)*vorz(il);
            }
            (ny<0)?(image2D[nz].negative[abs(ny)]=suml):(image2D[nz].negative[ny]=suml); //if ny<0 safe in negative part else (if ny>=0 safe in positive part)
        }
    }
}

void Multipole::doxyzimage(dataType grid)
{
    //TODO test if everything works fine with vector structure
    int max=63;
    image3D.resize(max);
    for(int nz=0;nz<=max;nz++)
    {
        dataType z=nz*grid;
        image3D[nz].expand(max);
        for(int ny=-max;ny<=max;ny++)
        {
            dataType y=ny*grid;
            //(ny<0)? (image3D[nz].negative[-ny].expand(max)):(image3D[nz].positive[ny].expand(max)); //expand only the enty ny that is used next
            SpecialVector<dataType> xtemp;
            xtemp.expand(max);
            for(int nx=-max;nx<=max;nx++)
            {
                dataType x=nx*grid;
                dataType r=sqrt(x*x+y*y+z*z);
                dataType alpha=acos(z/r);
                dataType beta=calcphi(x,y);
                dataType kr =k*r;
                dataType suml=0;
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


void Multipole::scaleimage(dataType grid)
{
    int max=63;
    assert(image3D.size()!=0);
    for(int iz=0;iz<=max;iz++)
    {
        for(int iy=-max;iy<max;iy++)
        {
            SpecialVector<dataType> xtemp;
            xtemp.expand(max);
            for(int ix=-max;ix<max;ix++)
            {
                dataType r=grid*sqrt(dataType(ix*ix+iy*iy+iz*iz));
                dataType temp= r*image3D[iz].getElement(iy).getElement(ix);
                xtemp.addElement(temp, ix);

            }
            image3D[iz].addElement(xtemp,iy);
        }
    }

}

void Multipole::smooth(dataType grid)
{
   /* assert(image3D.size()!=0);
    int box = 7; //smoothing box size half
    dataType width = 1; //gauss sigma width
    std::cout<<"pewforming gaussian smoothing\n";
    std::vector<SpecialVector<SpecialVector<dataType> > > >simage;
    simage=image3D;
    for(int nz=0;nz<=63;nz++)
    {
        for(int ny=0;ny<=63;ny++)
        {
            for(int nx=0;nx<=63;nx++)
            {
                dataType imageR=image3D[nz].getElement(ny).getElement(nx);
                for(iz=-box;iz<=box;iz++)
                {
                    for(iy=-box;iy<=box;iy++)
                    {
                        for(ix=-box;ix<=box;ix++)
                        {
                            dataType distsq=(ix+iy+iz)*(ix+iy+iz);
                            dataType gs=1/(width*2.5)*exp(-0.5*distsq/(width*width));
                            imageR+=sigma
                        }
                    }
                }
            }
        }
    }*--->/
}


/*-------------------Private--------------------------------------*/
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


/*dataType Multipole::calcth(dataType y, dataType z)
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

dataType Multipole::calcphi(dataType x, dataType y)
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
