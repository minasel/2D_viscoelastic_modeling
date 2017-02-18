! ********************************************************************** C
!                                                                        C
!  Program  :  getfai.                                                   C
!  Coded by :  Tong Xu                                                   C
!  Date     :  1995                                                      C
!  Language :  Fortran 77                                                C
!  Copyright:  Center for Lithospheric Studies                           C
!              The University of Texas at Dallas, 1995                   C
!                                                                        C
! ********************************************************************** C
! ======================================================================
      subroutine getfai_vsc(vp,vs,rho,qp,qs1,fai1,fai2,vpe,vse)
      implicit none
      include "commons.h"

      real vpe(-lv:nz+lv,-lv:nx+lv),vse(-lv:nz+lv,-lv:nx+lv)
      real vp(-lv:nz+lv,-lv:nx+lv),vs(-lv:nz+lv,-lv:nx+lv)
      real rho(-lv:nz+lv,-lv:nx+lv)
      real qp(-lv:nz+lv,-lv:nx+lv),qs1(-lv:nz+lv,-lv:nx+lv)
      real fai1(-lv:nz+lv,-lv:nx+lv,nq),fai2(-lv:nz+lv,-lv:nx+lv,nq)
!  Note:
! 
!      after this routine:
!
!    vp --->  lambda*dt/dx
!    vs --->  mu*dt/dx
!    
!
!  Arguments
      integer ix, iz, i, k,iq
      real sumfai1,sumfai2
!
      dtdx = dt/h
!

!
         !use qp to store relaxed compressional moduli
         !use qs1 to store relaxed shear moduli
         !use vp to store un-relaxed compressional moduli
         !use vs to store un-relaxed shear moduli
         do ix = -lv,nx+lv
         do iz = -lv,nz+lv
            qp(iz,ix) = rho(iz,ix)*vp(iz,ix)**2
            qs1(iz,ix) = rho(iz,ix)*vs(iz,ix)**2
            sumfai1 = 0.
            sumfai2 = 0.
            do iq=1, nq
               sumfai1 = sumfai1+fai1(iz,ix,iq)
               sumfai2 = sumfai2+fai2(iz,ix,iq)
            enddo
            vp(iz,ix) = qp(iz,ix)*(1-sumfai1)
            vs(iz,ix) = qs1(iz,ix)*(1-sumfai2)
         enddo
         enddo

         !fai1 (for dilatation),fai2 (for quasi-shear) 
         !are response functions from relq 
         !and complex modulus
         do iq = 1,nq
            do ix = -lv,nx+lv
            do iz = -lv,nz+lv
               fai1(iz,ix,iq) = qp(iz,ix)*fai1(iz,ix,iq)*relq(iq)
               fai2(iz,ix,iq) = qs1(iz,ix)*fai2(iz,ix,iq)*relq(iq)
            enddo
            enddo
         enddo

         !mu from unrelaxed shear moduli
         !now vs is mu*dt/dx
         do i=-lv,nx+lv
         do k=-lv,nz+lv
            vs(k,i) = dtdx*vs(k,i)
         end do
         end do

         !lambda from unrelaxed compressional moduli 
         !and shear modulus
         !now vp is lambda*dt/dx
         do i=-lv,nx+lv
         do k=-lv,nz+lv
            vp(k,i) = dtdx*vp(k,i)- 2*vs(k,i)
         end do
         end do
!
      !rarity from density
      !now rho is dt/(dx*rho)
      do i = -lv, nx+lv
      do k = -lv, nz+lv
        rho(k,i) = dtdx/rho(k,i)
      end do
      end do
!
      return
      end
!                                       
!=======================================================================
