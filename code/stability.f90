! ********************************************************************** C
!                                                                        C
!  Program  :  stability.f                                               C
!  Coded by :  Feng Deng                                                 C
!  Date     :  02/16/2007                                                C
!  Language :  Fortran 77                                                C
!  Copyright:  Center for Lithospheric Studies                           C
!              The University of Texas at Dallas, 2007                   C
!                                                                        C
! ********************************************************************** C
!=======================================================================
!
      subroutine stabil(vp, vs)
!
!----------------------------------------------------------------C
!   subroutine stabil                                            C
!   coded by   : Tong Xu                                         C
!   revised by : Feng Deng                                       C 
!                                                                C
!   This subroutine checks the stability for finite-differencing C
!   extrapolation                                                C
!                                                                C
!   *input*                                                      C
!      maxz      --  maximum index number of z                   C
!      maxx      --  maximum index number of x                   C
!      nz        --  number of grid points of z                  C
!      nx        --  number of grid points of x                  C
!      h         --  spatial sampling interval (dx=dz=h)         C
!      dt        --  time sampling interval                      C
!      vp        --  P wave velocity distribution (grid)         C
!      vs        --  S wave velocity distribution (grid)         C
!                                                                C
!  *output*                                                      C
!      if not stable, this subroutine will stop the whole        C
!      program                                                   C
!----------------------------------------------------------------C
!

      implicit none
      include "commons.h"
!  Arguments
      real  vp(-lv:nz+lv,-lv:nx+lv), vs(-lv:nz+lv,-lv:nx+lv)
!  Local Parameters
      integer i, k
      real  vmin, vmax, dtmax

! check accuracy

      vmin = 1.00e30

      do i = 1, nx
      do k = 1, nz
        vmin = amin1(vmin, vp(k,i))
      end do
      end do
      do i = 1, nx
      do k = 1, nz
        if(vs(k,i).gt.0.0) vmin = amin1(vmin,vs(k,i))
      end do
      end do

      vmin = sqrt(vmin)

! check stability condition

      vmax = 0.0

      do i = 1, nx
      do k = 1, nz
        vmax = amax1(vmax,vp(k,i))
      end do
      end do

! find stability limit

      dtmax = 0.70 * h / vmax / sqrt(2.0)                              
      if (dt .gt. dtmax) stop ' dt too large'


      return
      end
!
!=======================================================================
