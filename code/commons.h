      real, parameter :: pi=3.141592654
      real, parameter :: real0=0.d0
      integer, parameter :: mac_len=4
! lv: fixed for 8th-order finite-difference
      integer, parameter :: lv=4  
! Staggered finite-difference coefficients
      real, parameter, dimension(4)::&
          aa = (/1.23538628085693269, &
                -0.109160358903484037, &
                 0.246384765606012697E-1, &
                -0.660472512788868975E-2/)


!Model Parameters
        real  freq0
        real  h, dt
        real  relq(3)
        integer nsr,nk, nq,itsnap
        integer nx,ny,nz,nt
        integer knx,knz
        integer nzsr,nxsr
        integer it
        integer isource
	integer snap
	integer npml
        integer itshow
    	integer idpmin,idpmax   
	character*128 datapath
	real ak
	real dtdx



!
      common /idpmin/idpmin
      common /idpmax/idpmax
      common /dt/dt
      common /freq0/freq0
      common /datapath/datapath
      common /h/h
      common /itsnap/itsnap
      common /snap/snap
      common /npml/npml
      common /it/it
      common /isource/isource
      common /relq/relq
      common /itshow/itshow
      common /dtdx/dtdx
      common /nsr/nsr
      common /nq/nq
!
!
!
      common  /nx/nx
      common  /nz/nz
      common  /nt/nt
      common  /nk/nk
      common  /knx/knx
      common  /knz/knz
      common  /nzsr/nzsr
      common  /nxsr/nxsr
!



