! ********************************************************************** C
!                                                                        C
!  Program  :  getmod                                                    C
!  Coded by :  Wenlong Wang                                              C
!  Function :  Subroutine to read in parameter file                      C
!  Date     :  2015                                                      C
!  Language :  Fortran 90                                                C
!  Library  :  FFWT                                                      C
!  Copyright:  Center for Lithospheric Studies                           C
!              The University of Texas at Dallas, 2015                   C
!                                                                        C
! ********************************************************************** C
! ======================================================================!
      subroutine readparameter(inpt,zsrg1, xsrg1,src_int)
!
!
      implicit none
      include "commons.h"

! Arguments
      integer inpt
      integer zsrg1, xsrg1,src_int
! Local Parameter
      integer isr, iq, i, j
!
      !initial
!
      !read in model parameters
      call skiplines(inpt, 3)
      read(inpt, *) nz, nx, nt
      read(inpt, *) dt, h
      read(inpt, *) npml
!
      !read in visco-related parameters
      call skiplines(inpt, 3)
      read(inpt, *) nq 
      read(inpt, *) (relq(iq), iq=1, nq)
        do iq=1,nq
                relq(iq)=2*pi*relq(iq)
        enddo
!
      !read in miscellaneou parameters
      call skiplines(inpt, 3)
      read(inpt, *) snap
      read(inpt, *) itsnap
      read(inpt, *) itshow
!
      !read in source locations (z,x)
      call skiplines(inpt, 3)
      read(inpt, *) nsr
      read(inpt, *) zsrg1, xsrg1,src_int
      read(inpt, *) isource

!
      !read in dominant frequencys
      read(inpt, *) freq0

!
      !read in data path
      call skiplines(inpt, 3)
      read(inpt,'(a)') datapath
!
      end
