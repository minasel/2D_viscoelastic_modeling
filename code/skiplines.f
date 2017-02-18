C ********************************************************************** C
C                                                                        C
C  Program  :  skiplines.f                                               C
C  Coded by :  Feng Deng                                                 C
C  Date     :  03/16/2006                                                C
C  Language :  Fortran 77                                                C
C  Copyright:  Center for Lithospheric Studies                           C
C              The University of Texas at Dallas, 2007                   C
C                                                                        C
C ********************************************************************** C
         subroutine skiplines(lfilnam,n)
!==========================================================================
! This subroutine simply read file n lines
!--------------------------------------------------------------------------
! Author : Houzhu Zhang
!          University of Texas at Dallas, zhangh@utdallas.edu
!          9/1/2000
!==========================================================================
         implicit none

         integer  lfilnam,i,n

         do i=1,n
            read(lfilnam,*)
         enddo

         return
         end

