C ********************************************************************** C
C                                                                        C
C  Program  :  num2char.f                                                C
C  Coded by :  Feng Deng                                                 C
C  Date     :  02/16/2007                                                C
C  Language :  Fortran 77                                                C
C  Copyright:  Center for Lithospheric Studies                           C
C              The University of Texas at Dallas, 2007                   C
C                                                                        C
C ********************************************************************** C
C
        subroutine ficnum(i,c4i)
C
C----------------------------------------------------------C
C    subroutine ficnum                                     C
C                                                          C
C        change integer number to character (integer       C
C    be smalled than 9999)                                 C
C    *input*                                               C
C       i    : integer                                     C
C    *output*                                              C
C       c4i  : character                                   C
C                                                          C
C----------------------------------------------------------C
C
        implicit none
C
        integer i
        character *4 c4i
        character *3 c3i
        character *2 c2i
        character *1 c1i

        if ((i.lt.1).or.(i.gt.9999)) then
           write(*,*)'ficnum : number problem',i
           stop
        endif
C
        if (i.lt.10) then
           write(c1i,'(i1)') i
           c4i=c1i//'   '
        else if (i.lt.100) then
           write(c2i,'(i2)') i
           c4i=c2i//'  '
        else if (i.lt.1000) then
           write(c3i, '(i3)') i
           c4i=c3i//' '
        else
           write(c4i,'(i4)') i
        endif

        return
        end
C
c----------------------------------------------------------------
c////////////////////////////////////////////////////////////////
c----------------------------------------------------------------
c
         subroutine centra(name,imin,imax)
C
C----------------------------------------------------------C
C     subroutine centra                                    C
C                                                          C
C        to get the non-blank character in a string        C
C     *input*                                              C
C        name    :  string                                 C
C     *output*
C        imin    :  minimum index                          C
C        imax    :  maximum index                          C
C                                                          C
C----------------------------------------------------------C
         implicit none

         integer   i,imin,imax
         character *128 name

         do 1 i=1,128
            if (name(i:i).ne.' ') goto 11
1        continue

11       imin=i
         do 2 i=1,128
            if (name(129-i:129-i).ne.' ') goto 22
2        continue

22       imax=129-i

         return
         end
C
C=======================================================================
