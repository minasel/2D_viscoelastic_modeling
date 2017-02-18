C ********************************************************************* C
c
c  This subroutine:
c
c      Creates a source wavelet as a function of time steps.  The source
c  is defined as the derivative of a Gaussian function with respect to
c  time.
c
c      Returns array "seri", containing source wavelet.
c
c **********************************************************************
c
c     Coded by Robert Sun
c     Center for Lithospheric Studies
c     University of Texas at Dallas
c     January 1989
c     Version 1.0
c     Language: Fortran 77
c     Copyright: The University of Texas at Dallas
c
c     Modified by Feng Deng
c
c **********************************************************************
c
c     Modification by Xinfa Zhu, Nov 15, 2011
c     nprd is not easily understandable for users
c     change nprd to f0, dominant frequency of the output wavelet
c
c     nprd = 1.0/(1.4*f0*dt)
c
c     f0 in Hz, dt in second, 
c     nprd is the number of time points per half period.    
c     the output is a 90-degree-phase Ricker wavelet, anti-symmetric,
c       starts from time zero (it = 1);
c       negative at the first half period (it = 1, nprd);
c       positive at the second half period (it = nprd+1, nprd*2);
c       zero for other times (it = nprd*2+1, nt)
c       maximum amplitude is 1.0, minimum amplitude is -1.0
c     Looks like a negative sine function from 0 to 360 degrees.
c     The above relation is found and tested by Xinfa Zhu
c
c **********************************************************************
c
c  Variable :
c
c   nstep  : # of time steps.
c
c **********************************************************************
c
c  Array :
c
c   seri    : source as function of time.
c
c **********************************************************************
c
c     subroutine series(seri, minstep, maxstep, dt, nstep, nprd)
      subroutine series(seri, minstep, maxstep, dt, nstep, fdom)
C
      implicit none
C  Arguments
      integer nstep, nprd, minstep, maxstep
      real  dt, seri(minstep:maxstep)
C  Local Parameters
      integer it
      real  pi, amax, fdom, a, t,t0
c
c     added by Xinfa Zhu
      nprd = nint(1.0/(1.4*fdom*dt))
C
      pi=4.*atan(1.)
c assign 0 to maximum value.
      amax=0.
c choose source as derivative (with respect to time) of a gaussian.
C Option 1
c      do 815 it=1,nstep
c       seri(it)=-2.*(it-nprd+2)*exp(-(pi*dt*(it-nprd+2)/(nprd*dt))**2)
c         if(abs(seri(it)).gt.amax)amax=abs(seri(it))
c 815  continue

C Option 2
      a=pi*pi*fdom*fdom
      t0 = 1.20 /fdom 
      do 818 it=1,nstep
        t = real(it-1)*dt
        seri(it) = (1.e0 - 2.e0*a*(t-t0)**2)*exp(-a*(t-t0)**2)
         if(abs(seri(it)).gt.amax)amax=abs(seri(it))
 818  continue
c
c normalize source function.
      do 817 it=1,nstep
         seri(it)=seri(it)/amax
 817  continue
c
      return
      end
C
C ********************************************************************* C
