! ********************************************************************** C
!                                                                        C
!  Program  :  p_src.f                                                   C
!  Coded by :  Wenlong Wang                                              C
!  Date     :  03/22/2014                                                C
!  Language :  Fortran 90                                                C
!  Copyright:  Center for Lithospheric Studies                           C
!              The University of Texas at Dallas, 2014                   C
!                                                                        C
! ********************************************************************** C
!
        subroutine  p_src(vx,vz,nzsr1, nxsr1, amp)

!
!----------------------------------------------------------C
!    subroutine p_src                                      C
!                                                          C
!        Add pure P-wave source on particle velocities     C
!    *input*                                               C
!       vx,vz    : particle velocities                     C
!       nzsr1,nxsr1: source position
!       amp      : source amplitude
!                                                          C
!----------------------------------------------------------C
!
        implicit none
        include "commons.h"
        integer nzsr1, nxsr1 
!  Arguments
      real vx(-lv:nz+lv,-lv:nx+lv),vz(-lv:nz+lv,-lv:nx+lv)
      real amp

            vz(nzsr1,nxsr1)=vz(nzsr1,nxsr1)-aa(1)*amp
            vz(nzsr1+1,nxsr1)=vz(nzsr1+1,nxsr1)+aa(1)*amp
            vx(nzsr1+1,nxsr1)=vx(nzsr1+1,nxsr1)+aa(1)*amp
            vx(nzsr1+1,nxsr1-1)=vx(nzsr1+1,nxsr1-1)-aa(1)*amp

            vz(nzsr1-1,nxsr1)=vz(nzsr1-1,nxsr1)-aa(2)*amp
            vz(nzsr1+2,nxsr1)=vz(nzsr1+2,nxsr1)+aa(2)*amp
            vx(nzsr1+1,nxsr1+1)=vx(nzsr1+1,nxsr1+1)+aa(2)*amp
            vx(nzsr1+1,nxsr1-2)=vx(nzsr1+1,nxsr1-2)-aa(2)*amp
            vz(nzsr1-2,nxsr1)=vz(nzsr1-2,nxsr1)-aa(3)*amp
            vz(nzsr1+3,nxsr1)=vz(nzsr1+3,nxsr1)+aa(3)*amp
            vx(nzsr1+1,nxsr1+2)=vx(nzsr1+1,nxsr1+2)+aa(3)*amp
            vx(nzsr1+1,nxsr1-3)=vx(nzsr1+1,nxsr1-3)-aa(3)*amp
            vz(nzsr1-3,nxsr1)=vz(nzsr1-3,nxsr1)-aa(4)*amp
            vz(nzsr1+4,nxsr1)=vz(nzsr1+4,nxsr1)+aa(4)*amp
            vx(nzsr1+1,nxsr1+3)=vx(nzsr1+1,nxsr1+3)+aa(4)*amp
            vx(nzsr1+1,nxsr1-4)=vx(nzsr1+1,nxsr1-4)-aa(4)*amp

        return
        end
!
!----------------------------------------------------------------
!////////////////////////////////////////////////////////////////
!----------------------------------------------------------------
!
