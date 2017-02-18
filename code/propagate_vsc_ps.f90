! ********************************************************************** C
!                                                                        C
!  Program  :  propagate_vsc_ps.f90                                      C
!  Coded by :  Wenlong Wang                                              C
!  Date     :  03/23/2014                                                C
!  Language :  Fortran 90                                                C
!  Copyright:  Center for Lithospheric Studies                           C
!              The University of Texas at Dallas, 2014                   C
!                                                                        C
! ********************************************************************** C
! ===================================================================
!
      subroutine eqmot_vis_ps(vp,vs,rho,vx,vz,sxx,szz,szx,&
                       vpx,vpz,vsx,vsz,spp,&
                       memory_dsigmaxx_dx,memory_dsigmazz_dz,&
                       memory_dsigmazx_dx,memory_dsigmazx_dz,&
                       a_z,b_z,k_z,a_x,b_x,k_x,&
                   a_z_half, b_z_half,K_z_half,a_x_half,b_x_half,K_x_half)


!
! -------------------------------------------------------------------
! equation of motion
! compute particle velocities and displacement from stresses
! Copyright : The University of Texas at Dallas, 2014.
! -------------------------------------------------------------------
!
      implicit none
      include "commons.h"
!
!  Arguments
      real vp(-lv:nz+lv,-lv:nx+lv),vs(-lv:nz+lv,-lv:nx+lv)
      real rho(-lv:nz+lv,-lv:nx+lv)
      real vx(-lv:nz+lv,-lv:nx+lv),vz(-lv:nz+lv,-lv:nx+lv)
      real vsx(-lv:nz+lv,-lv:nx+lv),vsz(-lv:nz+lv,-lv:nx+lv)
      real vpx(-lv:nz+lv,-lv:nx+lv),vpz(-lv:nz+lv,-lv:nx+lv)
      real sxx(-lv:nz+lv,-lv:nx+lv),szz(-lv:nz+lv,-lv:nx+lv)
      real szx(-lv:nz+lv,-lv:nx+lv),spp(-lv:nz+lv,-lv:nx+lv)
      real memory_dsigmaxx_dx(-lv:nz+lv,-lv:nx+lv),memory_dsigmazz_dz(-lv:nz+lv,-lv:nx+lv)
      real memory_dsigmazx_dx(-lv:nz+lv,-lv:nx+lv),memory_dsigmazx_dz(-lv:nz+lv,-lv:nx+lv)
      real a_z(1:nz),b_z(1:nz),k_z(1:nz)
      real a_x(1:nx),b_x(1:nx),k_x(1:nx)
      real a_z_half(1:nz),b_z_half(1:nz),k_z_half(1:nz)
      real a_x_half(1:nx),b_x_half(1:nx),k_x_half(1:nx)
!  Local Parameters
      integer i, k
      real value_dsigmaxx_dx, value_dsigmazx_dz
      real value_dsigmazz_dz, value_dsigmazx_dx
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,value_dsigmaxx_dx),&
!$OMP& PRIVATE(value_dsigmazx_dz)

!$OMP DO SCHEDULE(DYNAMIC,1)
      do i = 1, nx
      do k = 1, nz
        value_dsigmaxx_dx=aa(1)*(sxx(k,i+1)-sxx(k,i  ))&
                         +aa(2)*(sxx(k,i+2)-sxx(k,i-1))&
                         +aa(3)*(sxx(k,i+3)-sxx(k,i-2))&
                         +aa(4)*(sxx(k,i+4)-sxx(k,i-3))
        value_dsigmazx_dz=aa(1)*(szx(k  ,i)-szx(k-1,i))&
                         +aa(2)*(szx(k+1,i)-szx(k-2,i))&
                         +aa(3)*(szx(k+2,i)-szx(k-3,i))&
                         +aa(4)*(szx(k+3,i)-szx(k-4,i))
      memory_dsigmaxx_dx(k,i) = b_x_half(i) * memory_dsigmaxx_dx(k,i) +&
                               a_x_half(i) * value_dsigmaxx_dx
      memory_dsigmazx_dz(k,i) = b_z(k) * memory_dsigmazx_dz(k,i) +&
                               a_z(k) * value_dsigmazx_dz
      value_dsigmaxx_dx = value_dsigmaxx_dx / K_x_half(i) +&
                          memory_dsigmaxx_dx(k,i)
      value_dsigmazx_dz = value_dsigmazx_dz / K_z(k) +&
                           memory_dsigmazx_dz(k,i)


        vx(k,i) = vx(k,i)&
         + (value_dsigmaxx_dx+value_dsigmazx_dz)&
         * (rho(k,i)+rho(k,i+1)+rho(k+1,i)+rho(k+1,i+1))*0.25

        vpx(k,i) = vpx(k,i)&
         + (aa(1)*(spp(k,i+1)-spp(k,i  ))&
           +aa(2)*(spp(k,i+2)-spp(k,i-1))&
           +aa(3)*(spp(k,i+3)-spp(k,i-2))&
           +aa(4)*(spp(k,i+4)-spp(k,i-3)))&
           *(rho(k,i)+rho(k,i+1)+rho(k+1,i)+rho(k+1,i+1))*0.25
        vsx(k,i) = vx(k,i)-vpx(k,i)

      end do
      end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,value_dsigmazx_dx,value_dsigmazz_dz)

!$OMP DO SCHEDULE(DYNAMIC,1)

      do i = 1, nx
      do k = 1, nz
        value_dsigmazx_dx=aa(1)*(szx(k,i  )-szx(k,i-1))&
                         +aa(2)*(szx(k,i+1)-szx(k,i-2))&
                         +aa(3)*(szx(k,i+2)-szx(k,i-3))&
                         +aa(4)*(szx(k,i+3)-szx(k,i-4))
        value_dsigmazz_dz=aa(1)*(szz(k+1,i)-szz(k  ,i))&
                         +aa(2)*(szz(k+2,i)-szz(k-1,i))&
                         +aa(3)*(szz(k+3,i)-szz(k-2,i))&
                         +aa(4)*(szz(k+4,i)-szz(k-3,i))
      memory_dsigmazx_dx(k,i) = b_x(i) * memory_dsigmazx_dx(k,i) +&
                               a_x(i) * value_dsigmazx_dx
      memory_dsigmazz_dz(k,i) = b_z_half(k) * memory_dsigmazz_dz(k,i) +&
                               a_z_half(k) * value_dsigmazz_dz
      value_dsigmazx_dx = value_dsigmazx_dx / K_x(i) +&
                         memory_dsigmazx_dx(k,i)
      value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(k) +&
                         memory_dsigmazz_dz(k,i)

        vz(k,i) = vz(k,i)&
         + (value_dsigmazx_dx+value_dsigmazz_dz)&
          *(rho(k,i)+rho(k,i+1)+rho(k+1,i)+rho(k+1,i+1))*0.25

        vpz(k,i) = vpz(k,i)&
          + (aa(1)*(spp(k+1,i)-spp(k  ,i))&
            +aa(2)*(spp(k+2,i)-spp(k-1,i))&
            +aa(3)*(spp(k+3,i)-spp(k-2,i))&
            +aa(4)*(spp(k+4,i)-spp(k-3,i)))&
          *(rho(k,i)+rho(k,i+1)+rho(k+1,i)+rho(k+1,i+1))*0.25

        vsz(k,i) = vz(k,i)-vpz(k,i)

      end do
      end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL
!
      return
      end

!-----------------------------------------------------------------
!/////////////////////////////////////////////////////////////////
!-----------------------------------------------------------------
!
      subroutine hooke_vis_ps(vp,vs,rho,vx,vz,sxx,szz,szx,spp,&
                       qp,qs1,qs2,divergv,&
                       memory_dvx_dx,memory_dvz_dz,&
                       memory_dvz_dx,memory_dvx_dz,&
                       a_z,b_z,k_z,a_x,b_x,k_x,&
                       a_z_half, b_z_half,K_z_half,a_x_half,b_x_half,K_x_half)

 
! -------------------------------------------------------------------
! generalized hooke's law
! compute stresses from particle velocities, displacement and memory
! variables
! Copyright : The University of Texas at Dallas, 2014.
! -------------------------------------------------------------------
      implicit none
      include "commons.h"
!
!  Arguments
      real vp(-lv:nz+lv,-lv:nx+lv),vs(-lv:nz+lv,-lv:nx+lv)
      real rho(-lv:nz+lv,-lv:nx+lv)
      real vx(-lv:nz+lv,-lv:nx+lv),vz(-lv:nz+lv,-lv:nx+lv)
      real sxx(-lv:nz+lv,-lv:nx+lv),szz(-lv:nz+lv,-lv:nx+lv)
      real szx(-lv:nz+lv,-lv:nx+lv),spp(-lv:nz+lv,-lv:nx+lv)
      real qs1(-lv:nz+lv,-lv:nx+lv),qp(-lv:nz+lv,-lv:nx+lv)
      real qs2(-lv:nz+lv,-lv:nx+lv),divergv(-lv:nz+lv,-lv:nx+lv)
      real memory_dvx_dx(-lv:nz+lv,-lv:nx+lv),memory_dvz_dz(-lv:nz+lv,-lv:nx+lv)
      real memory_dvz_dx(-lv:nz+lv,-lv:nx+lv),memory_dvx_dz(-lv:nz+lv,-lv:nx+lv)
      real a_z(1:nz),b_z(1:nz),k_z(1:nz)
      real a_x(1:nx),b_x(1:nx),k_x(1:nx)
      real a_z_half(1:nz),b_z_half(1:nz),k_z_half(1:nz)
      real a_x_half(1:nx),b_x_half(1:nx),k_x_half(1:nx)

!  Local Parameters
      integer i, k
      real*4 value_dvx_dx, value_dvz_dz
      real*4 value_dvz_dx, value_dvx_dz

! vp = lambda*dt/h
! vs =     mu*dt/h


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,value_dvx_dx,value_dvz_dz)

!$OMP DO SCHEDULE(DYNAMIC,1)
      do i = 1, nx
      do k = 1, nz
        value_dvx_dx=aa(1)*(vx(k,i)-vx(k,i-1))&
                    +aa(2)*(vx(k,i+1)-vx(k,i-2))&
                    +aa(3)*(vx(k,i+2)-vx(k,i-3))&
                    +aa(4)*(vx(k,i+3)-vx(k,i-4))
        value_dvz_dz=aa(1)*(vz(k,i)-vz(k-1,i))&
                    +aa(2)*(vz(k+1,i)-vz(k-2,i))&
                    +aa(3)*(vz(k+2,i)-vz(k-3,i))&
                    +aa(4)*(vz(k+3,i)-vz(k-4,i))
        memory_dvx_dx(k,i) = b_x(i) * memory_dvx_dx(k,i) + a_x(i) *&
                             value_dvx_dx
        memory_dvz_dz(k,i) = b_z(k) * memory_dvz_dz(k,i) + a_z(k) *&
                             value_dvz_dz
        value_dvx_dx = value_dvx_dx / K_x(i) + memory_dvx_dx(k,i)
        value_dvz_dz = value_dvz_dz / K_z(k) + memory_dvz_dz(k,i)


        sxx(k,i) = sxx(k,i)&
           +(value_dvx_dx+value_dvz_dz)&
           *(vp(k,i)+2.*vs(k,i)) - value_dvz_dz*2*vs(k,i)&
           +(qp(k,i)+qs1(k,i))*dt

        szz(k,i) = szz(k,i)&
           +(value_dvx_dx+value_dvz_dz)&
           *(vp(k,i)+2.*vs(k,i)) - value_dvx_dx*2*vs(k,i)&
           +(qp(k,i)+qs2(k,i))*dt

        spp(k,i) = spp(k,i)&
           + (aa(1)*(vz(k  ,i)-vz(k-1,i))+aa(2)*(vz(k+1,i)-vz(k-2,i))&
             +aa(3)*(vz(k+2,i)-vz(k-3,i))+aa(4)*(vz(k+3,i)-vz(k-4,i))&
             +aa(1)*(vx(k,i  )-vx(k,i-1))+aa(2)*(vx(k,i+1)-vx(k,i-2))&
             +aa(3)*(vx(k,i+2)-vx(k,i-3))+aa(4)*(vx(k,i+3)-vx(k,i-4)))&
            *(vp(k,i)+2.0*vs(k,i))+qp(k,i)*dt

      end do
      end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,value_dvz_dx,value_dvx_dz)

!$OMP DO SCHEDULE(DYNAMIC,1)
      do i = 1, nx
      do k = 1, nz
        value_dvz_dx=aa(1)*(vz(k,i+1)-vz(k,i))&
                    +aa(2)*(vz(k,i+2)-vz(k,i-1))&
                    +aa(3)*(vz(k,i+3)-vz(k,i-2))&
                    +aa(4)*(vz(k,i+4)-vz(k,i-3))
        value_dvx_dz=aa(1)*(vx(k+1,i)-vx(k,i))&
                    +aa(2)*(vx(k+2,i)-vx(k-1,i))&
                    +aa(3)*(vx(k+3,i)-vx(k-2,i))&
                    +aa(4)*(vx(k+4,i)-vx(k-3,i))
     memory_dvz_dx(k,i)=b_x_half(i)*memory_dvz_dx(k,i)+a_x_half(i)*&
                          value_dvz_dx
     memory_dvx_dz(k,i)=b_z_half(k)*memory_dvx_dz(k,i)+a_z_half(k)*&
                          value_dvx_dz

        value_dvz_dx = value_dvz_dx / K_x_half(i) + memory_dvz_dx(k,i)
        value_dvx_dz = value_dvx_dz / K_z_half(k) + memory_dvx_dz(k,i)

        szx(k,i) = szx(k,i)&
          +(value_dvz_dx+value_dvx_dz)*vs(k,i)&
          +divergv(k,i)*dt
      end do
      end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL
!
      return
      end
!
!-----------------------------------------------------------------
!/////////////////////////////////////////////////////////////////
!-----------------------------------------------------------------
!
!
      subroutine memory_ps(vx,vz,memory_dvx_dx,memory_dvz_dz,&
                   memory_dvx_dz,memory_dvz_dx,k_x,k_z,k_x_half,k_z_half,&
                   memory_qs1,memory_qs2,fai1,fai2,qp,qs1,qs2,divergv,&
                   memory_xz,memory_qp)
!
! -------------------------------------------------------------------
! use Carcione's memory variables to get viscoelastic results
! Copyright : The University of Texas at Dallas, 2014.
! -------------------------------------------------------------------
!
      implicit none
      include "commons.h"
!
!  Arguments
      real vx(-lv:nz+lv,-lv:nx+lv),vz(-lv:nz+lv,-lv:nx+lv)
      real memory_dvz_dz(-lv:nz+lv,-lv:nx+lv),memory_dvx_dx(-lv:nz+lv,-lv:nx+lv)
      real memory_dvx_dz(-lv:nz+lv,-lv:nx+lv),memory_dvz_dx(-lv:nz+lv,-lv:nx+lv)
      real qs1(-lv:nz+lv,-lv:nx+lv),qs2(-lv:nz+lv,-lv:nx+lv),qp(-lv:nz+lv,-lv:nx+lv)
      real divergv(-lv:nz+lv,-lv:nx+lv)
      real memory_qp(-lv:nz+lv,-lv:nx+lv,nq)
      real memory_qs1(-lv:nz+lv,-lv:nx+lv,nq),fai2(-lv:nz+lv,-lv:nx+lv,nq)
      real memory_qs2(-lv:nz+lv,-lv:nx+lv,nq),fai1(-lv:nz+lv,-lv:nx+lv,nq)
      real memory_xz(-lv:nz+lv,-lv:nx+lv,nq)
      real k_z(1:nz),k_x(1:nz)
      real k_z_half(1:nz),k_x_half(1:nz)
!  Local Parameters
      integer l, i, k
      real*4  value_dvx_dx, value_dvz_dz, value_dvx_dz, value_dvz_dx
      real*4  dmemory_qp, dmemory_qs1, dmemory_qs2, dmemory_xz


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,l,value_dvx_dx,value_dvz_dz,dmemory_qp)

!$OMP DO SCHEDULE(DYNAMIC,1) COLLAPSE (2)

      do l = 1, nq
      do i = 1, nx
      do k = 1, nz
        value_dvx_dx = aa(1)*(vx(k,i  )-vx(k,i-1))&
                      +aa(2)*(vx(k,i+1)-vx(k,i-2))&
                      +aa(3)*(vx(k,i+2)-vx(k,i-3))&
                      +aa(4)*(vx(k,i+3)-vx(k,i-4))
        value_dvz_dz = aa(1)*(vz(k,i  )-vz(k-1,i))&
                      +aa(2)*(vz(k+1,i)-vz(k-2,i))&
                      +aa(3)*(vz(k+2,i)-vz(k-3,i))&
                      +aa(4)*(vz(k+3,i)-vz(k-4,i))

        value_dvx_dx = value_dvx_dx / K_x(i) + memory_dvx_dx(k,i)
        value_dvz_dz = value_dvz_dz / K_z(k) + memory_dvz_dz(k,i)

        dmemory_qp = (value_dvx_dx + value_dvz_dz)/h*fai1(k,i,l)
        memory_qp(k,i,l)=memory_qp(k,i,l)+dt*(dmemory_qp&
                       -memory_qp(k,i,l)*relq(l))


        value_dvx_dx = aa(1)*(vx(k,i  )-vx(k,i-1))&
                      +aa(2)*(vx(k,i+1)-vx(k,i-2))&
                      +aa(3)*(vx(k,i+2)-vx(k,i-3))&
                      +aa(4)*(vx(k,i+3)-vx(k,i-4))
        value_dvz_dz = aa(1)*(vz(k,i  )-vz(k-1,i))&
                      +aa(2)*(vz(k+1,i)-vz(k-2,i))&
                      +aa(3)*(vz(k+2,i)-vz(k-3,i))&
                      +aa(4)*(vz(k+3,i)-vz(k-4,i))

        value_dvx_dx = value_dvx_dx / K_x(i) + memory_dvx_dx(k,i)
        value_dvz_dz = value_dvz_dz / K_z(k) + memory_dvz_dz(k,i)

        dmemory_qs1 = (value_dvz_dz)/h*fai2(k,i,l)*2
        dmemory_qs2 = (value_dvx_dx)/h*fai2(k,i,l)*2

        memory_qs1(k,i,l)=memory_qs1(k,i,l)-dt*(dmemory_qs1 &
                         +memory_qs1(k,i,l)*relq(l))
        memory_qs2(k,i,l)=memory_qs2(k,i,l)-dt*(dmemory_qs2 &
                         +memory_qs2(k,i,l)*relq(l))

        value_dvx_dz = aa(1)*(vx(k+1,i)-vx(k  ,i))&
                      +aa(2)*(vx(k+2,i)-vx(k-1,i))&
                      +aa(3)*(vx(k+3,i)-vx(k-2,i))&
                      +aa(4)*(vx(k+4,i)-vx(k-3,i))
        value_dvz_dx = aa(1)*(vz(k,i+1)-vz(k,i  ))&
                      +aa(2)*(vz(k,i+2)-vz(k,i-1))&
                      +aa(3)*(vz(k,i+3)-vz(k,i-2))&
                      +aa(4)*(vz(k,i+4)-vz(k,i-3))

        value_dvz_dx = value_dvz_dx / K_x_half(i) + memory_dvz_dx(k,i)
        value_dvx_dz = value_dvx_dz / K_z_half(k) + memory_dvx_dz(k,i)

        dmemory_xz =  (value_dvz_dx+value_dvx_dz)/h*fai2(k,i,l)

        memory_xz(k,i,l) = memory_xz(k,i,l)+dt*(dmemory_xz &
                         -memory_xz(k,i,l)*relq(l))

!                                       
      end do
      end do
      end do

!$OMP END DO

!$OMP END PARALLEL
      do i = 1, nx
      do k = 1, nz
        qp(k,i) = 0.
        qs1(k,i) = 0.
        qs2(k,i) = 0.
        divergv(k,i)=0.
         do l=1, nq
            qp(k,i)  = qp(k,i)  + memory_qp(k,i,l)
            qs1(k,i) = qs1(k,i) + memory_qs1(k,i,l)
            qs2(k,i) = qs2(k,i) + memory_qs2(k,i,l)
            divergv(k,i)  = divergv(k,i)+memory_xz(k,i,l)
         enddo
      end do
      end do


      return
      end
!
!-----------------------------------------------------------------
!/////////////////////////////////////////////////////////////////
!-----------------------------------------------------------------

