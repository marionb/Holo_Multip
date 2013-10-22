#include <iostream>
#include <cmath>

#include "Multipole.h"
#include "Data.h"

/*
void multi()
 thetst=messg(index(1),2)-messg(index(2),2)
      bnorm=0.0
      dnorm=0.0
      do i=1,nt
         theta=messg(index(i),2)
         do m=0,lmax,isym
            rm=m
            rint1=0.0
            rint2=0.0
            do j=index(i),index(i+1)-1
               phi=messg(j,3)
               rint1=rint1+messg(j,1)*cosd(rm*phi)
               rint2=rint2+messg(j,1)*sind(rm*phi)
            end do
            bm1(i,m)=rint1*sind(theta)*dphi(i)*thetst
            bm2(i,m)=rint2*sind(theta)*dphi(i)*thetst
         end do
         bnorm=bnorm+bm1(i,0)
      end do

c     in our normalization bnorm is very near to zero, and
c     it would be dangerous to divide by bnorm
      if((inorm.eq.1).or.(bnorm.lt.0.001))bnorm=1.0

      do l=0,lmax,2
         il=(l/2)
         do m=0,l,isym
            aint=0.0
            bint=0.0
            do i=1,nt
               theta=messg(index(i),2)
               var=cosd(theta)
               algndr=plgndr(l,m,var)
               aint=aint+bm1(i,m)*algndr
               bint=bint-bm2(i,m)*algndr
            end do
            alm1(il,m)=aint/bnorm
            alm2(il,m)=bint/bnorm
         end do
      end do
      return
      end
}*/

int main()
{
    double grid=0.09;
    Multipole *file=new Multipole("oldinp.itp",50,1,3.73);
    std::cout<<"created instance\n";
    file->readData();
    file->multpl();
    //ile->printAlm();
    file->expans();

    file->writeData("testOut.tex");
    //file->doyzimage(grid,0);

    //std::cout<<"apofact is "<<apofct(file->messg[0][1],file->messg[0][1]);

    return 0;
}








