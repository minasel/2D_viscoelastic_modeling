!
      subroutine getmod(vpfile, vsfile, rhofile, vp, vs, rho, fai1, fai2)
!
!----------------------------------------------------------------C
!                                                                C
!    subroutine getmod                                           C
!                                                                C
!       to read in vp, vs, and density response functions        C
!                                                                C
!----------------------------------------------------------------C
!
      implicit none
      include "commons.h"
!  Arguments
      real fai1(-lv:nz+lv,-lv:nx+lv,nq),fai2(-lv:nz+lv,-lv:nx+lv,nq)
      real vp(-lv:nz+lv,-lv:nx+lv),vs(-lv:nz+lv,-lv:nx+lv)
      real rho(-lv:nz+lv,-lv:nx+lv)

!  Parameters
      character*128 vpfile, vsfile, rhofile
!  Local Parameters
      integer iz, ix, inpt, imin, imax, l, i, j, m
      character*128 filename
      character*4 cm
!
      !initial
      inpt = 61
!
      !read in p-velocity model
      call centra(vpfile,imin,imax)
      open(inpt, file=vpfile(imin:imax), access='direct',&
          form='unformatted', recl=mac_len*nz)
      do ix=1,nx
         read(inpt,rec=ix)(vp(iz,ix),iz=1,nz)
      enddo
      do iz=-lv,nz+lv
        do ix=-lv,0
                vp(iz,ix)=vp(iz,1)
        enddo
        do ix=nx,nx+lv
                vp(iz,ix)=vp(iz,nx)
         enddo
      enddo
      do ix=-lv,nx+lv
        do iz=-lv,0
                vp(iz,ix)=vp(1,ix)
        enddo
        do iz=nz,nz+lv
                vp(iz,ix)=vp(nz,ix)
        enddo
      enddo
      close(inpt)
!
      !read in s-velocity model
      call centra(vsfile,imin,imax)
      open(inpt, file=vsfile(imin:imax), access='direct',&
          form='unformatted', recl=mac_len*nz)
      do ix=1,nx
         read(inpt,rec=ix)(vs(iz,ix),iz=1,nz)
      enddo
      do iz=-lv,nz+lv
        do ix=-lv,0
                vs(iz,ix)=vs(iz,1)
        enddo
        do ix=nx,nx+lv
                vs(iz,ix)=vs(iz,nx)
        enddo
      enddo
      do ix=-lv,nx+lv
        do iz=-lv,0
                vs(iz,ix)=vs(1,ix)
        enddo
        do iz=nz,nz+lv
                vs(iz,ix)=vs(nz,ix)
        enddo
      enddo
      close(inpt)
!
      !read in density
      call centra(rhofile,imin,imax)
      open(inpt, file=rhofile(imin:imax), access='direct',&
          form='unformatted', recl=mac_len*nz)
      do ix=1,nx
         read(inpt,rec=ix)(rho(iz,ix),iz=1,nz)
      enddo
      do iz=-lv,nz+lv
        do ix=-lv,0
                rho(iz,ix)=rho(iz,1)
        enddo
        do ix=nx,nx+lv
                rho(iz,ix)=rho(iz,nx)
        enddo
      enddo
      do ix=-lv,nx+lv
        do iz=-lv,0
                rho(iz,ix)=rho(1,ix)
        enddo
        do iz=nz,nz+lv
                rho(iz,ix)=rho(nz,ix)
        enddo
      enddo
      close(inpt)
!
      ! read in relaxation times
            do m=1, nq
!
               call ficnum(m, cm)
               filename = datapath(idpmin:idpmax)//"tau_model_1"//cm
               call centra(filename, imin, imax)
               open(88, file=filename(imin:imax), access='direct',&
                   form='unformatted', recl=mac_len*nz)
               do ix=1, nx
                  read(88, rec=ix)(fai1(iz,ix,m), iz=1,nz)
               enddo
              do iz=-lv,nz+lv
                do ix=-lv,0
                        fai1(iz,ix,m)=fai1(iz,1,m)
                enddo
                do ix=nx,nx+lv
                        fai1(iz,ix,m)=fai1(iz,nx,m)
                enddo
              enddo
              do ix=-lv,nx+lv
                do iz=-lv,0
                        fai1(iz,ix,m)=fai1(1,ix,m)
                enddo
                do iz=nz,nz+lv
                        fai1(iz,ix,m)=fai1(nz,ix,m)
                enddo
              enddo
               close(88)
              do iz=-lv,nz+lv
              do ix=-lv,nx+lv
                        fai1(iz,ix,m)=1.-fai1(iz,ix,m)*relq(m)
              enddo
              enddo
!
!
               filename = datapath(idpmin:idpmax)//"tau_model_2"//cm
               call centra(filename, imin, imax)
               open(88, file=filename(imin:imax), access='direct',&
                   form='unformatted', recl=mac_len*nz)
               do ix=1, nx
                  read(88, rec=ix)(fai2(iz,ix,m), iz=1,nz)
               enddo
              do iz=-lv,nz+lv
                do ix=-lv,0
                        fai2(iz,ix,m)=fai2(iz,1,m)
                enddo
                do ix=nx,nx+lv
                        fai2(iz,ix,m)=fai2(iz,nx,m)
                enddo
              enddo
              do ix=-lv,nx+lv
                do iz=-lv,0
                        fai2(iz,ix,m)=fai2(1,ix,m)
                enddo
                do iz=nz,nz+lv
                        fai2(iz,ix,m)=fai2(nz,ix,m)
                enddo
              enddo
               close(88)
              do iz=-lv,nz+lv
              do ix=-lv,nx+lv
                        fai2(iz,ix,m)=1.-fai2(iz,ix,m)*relq(m)
              enddo
              enddo
!
            enddo
!
      return
      end
!
