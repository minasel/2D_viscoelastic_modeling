! Random change
! ********************************************************************** !
!                                                                        !
!  Program  :  velvob2d.f                                                !
!  CODE  BY :  Wenlong Wang                                              !
!  Date     :  03/18/2016                                                !
!  Language :  Fortran 90, with MPI                                      !
!  Copyright:  Center for Lithospheric Studies                           !
!              The University of Texas at Dallas, 2016                   !
!                                                                        !
! ********************************************************************** !
!
!    This is a parallelized (by MPI) over shots finite-difference 
!  viscoelastic modeling + PS decomposition program.
!
!      1) Generates the vertical and horizontal particle velocity seismogram 
!  which is the seismic response to the input model by viscoelastic
!  finite differences. these seismic data are used as the observed test 
!  data in the migration.
!      
!      2) Generates pure P and mixed P + S sources synthetic data.
!
!      3) Generates decomposed P- and S-wave snapshots.
!
!
! **********************************************************************
!
!  *Input*
!     
!    comp.par   : parameter file
!
!    vp_model    : P wave velocity model
!
!    vs_model    : S wave velocity model
!
!    rho_model   : density model
!
!    tau_model_ij    : stress stain relaxation time
!
!
!  *Output*
!   
!    f110#      : observed horizontal component of the particle velocity 
!                 of sesimic data
!
!    f120#      : observed vertical component of the particle velocity
!                 of sesimic data

!    s010#      : horizontal component of the complete snapshot
!
!    s020#      : vertical component of the complete snapshot

!    s110#      : horizontal component of the decomposed P-wave snapshot
!
!    s120#      : vertical component of the decomposed P-wave snapshot
!
!    s210#      : horizontal component of the decomposed S-wave snapshot
!
!    s220#      : vertical component of the decomposed S-wave snapshot
!
!    source.bin : source wavelet
!
!
!
! **********************************************************************
!
!  get space:
!
      use omp_lib
      implicit none
      include "commons.h"
      include "commons_alc.h"
      include "mpif.h"
! Local Parameters
      integer isr,  iq
      integer i, j, k, l, imin, imax, ix,iz
      integer ifoh, ifov, ivalid
      integer zsrg1,xsrg1,src_int, rec_z
      character*4  cisr
      integer isnap
      real amp
! Input File Names
      character*128 vpfile, vsfile, rhofile,filename
! MPI Parameters
      integer ierr, myid, numprcs, namelen
      double precision starttime, endtime
      character*(MPI_MAX_PROCESSOR_NAME) processorname
!
      !MPI Initialization
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprcs, ierr)
      call MPI_Get_processor_name(processorname, namelen, ierr)
!
      !output some node information
      write(0,301) myid, processorname(1:namelen)
      if(myid .eq. 0) then
         write(0,*) "the number of processes is:", numprcs
        
      endif
!
      !synchronize at this point and start to time the process
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      starttime = MPI_WTIME()
      write(0,*) "my id is ", myid
!
      !define I/O
      ifoh = 21
      ifov = 22
!
      !read input parameters.
      open(11, file='comp.par', status='old')
      call readparameter(11,zsrg1, xsrg1,src_int)
      close(11)

!       allocate source positions
      allocate(zsrg(nsr),xsrg(nsr))
      do i=1,nsr
              zsrg(i)=zsrg1
              xsrg(i)=xsrg1+src_int*(i-1)
      enddo

!! allocate memory for wavefield variables
!
!  Model related memory
        allocate(seri(nt))! nt is the length of the source wavelet, 
!       and must be the same as the number of time steps computed
! lv is fixed to be 4 for 8th-order FD
        allocate(vp(-lv:nz+lv,-lv:nx+lv),vs(-lv:nz+lv,-lv:nx+lv))
        allocate(rho(-lv:nz+lv,-lv:nx+lv))
        allocate(vpe(-lv:nz+lv,-lv:nx+lv),vse(-lv:nz+lv,-lv:nx+lv))
        allocate(qs2(-lv:nz+lv,-lv:nx+lv),qs1(-lv:nz+lv,-lv:nx+lv))
        allocate(qp(-lv:nz+lv,-lv:nx+lv),divergv(-lv:nz+lv,-lv:nx+lv))

        allocate(sxx(-lv:nz+lv,-lv:nx+lv),szz(-lv:nz+lv,-lv:nx+lv))
        allocate(szx(-lv:nz+lv,-lv:nx+lv),spp(-lv:nz+lv,-lv:nx+lv))

        allocate(vx(-lv:nz+lv,-lv:nx+lv),vz(-lv:nz+lv,-lv:nx+lv))
        allocate(vpx(-lv:nz+lv,-lv:nx+lv),vpz(-lv:nz+lv,-lv:nx+lv))
        allocate(vsx(-lv:nz+lv,-lv:nx+lv),vsz(-lv:nz+lv,-lv:nx+lv))
        allocate(fai1(-lv:nz+lv,-lv:nx+lv,nq),fai2(-lv:nz+lv,-lv:nx+lv,nq))
        allocate(memory_qp(-lv:nz+lv,-lv:nx+lv,nq),memory_qs1(-lv:nz+lv,-lv:nx+lv,nq))
        allocate(memory_qs2(-lv:nz+lv,-lv:nx+lv,nq),memory_xz(-lv:nz+lv,-lv:nx+lv,nq))


        !   PML vaiable allocation
        allocate(memory_dvx_dx(-lv:nz+lv,-lv:nx+lv),memory_dvz_dz(-lv:nz+lv,-lv:nx+lv))
        allocate(memory_dvx_dz(-lv:nz+lv,-lv:nx+lv),memory_dvz_dx(-lv:nz+lv,-lv:nx+lv))
        allocate(memory_dsigmaxx_dx(-lv:nz+lv,-lv:nx+lv),memory_dsigmazz_dz(-lv:nz+lv,-lv:nx+lv))
        allocate(memory_dsigmazx_dx(-lv:nz+lv,-lv:nx+lv),memory_dsigmazx_dz(-lv:nz+lv,-lv:nx+lv))
        allocate(a_z(1:nz),b_z(1:nz),k_z(1:nz))
        allocate(a_x(1:nx),b_x(1:nx),k_x(1:nx))
        allocate(a_z_half(1:nz),b_z_half(1:nz),k_z_half(1:nz))
        allocate(a_x_half(1:nx),b_x_half(1:nx),k_x_half(1:nx))
!
!      initialize variables
!
        seri=real0
        vp=real0
        vs=real0
        rho=real0
        vpe=real0
        vse=real0
        fai1=real0
        fai2=real0
        a_z=real0
        b_z=real0
        k_z=real0
        a_x=real0
        b_x=real0
        k_x=real0
        a_z_half=real0
        b_z_half=real0
        k_z_half=real0
        a_x_half=real0
        b_x_half=real0
        k_x_half=real0
!
! Read in velocity models
!
      call centra(datapath, idpmin, idpmax)
      vpfile = datapath(idpmin:idpmax)//"vp_model"
      vsfile = datapath(idpmin:idpmax)//"vs_model"
      rhofile = datapath(idpmin:idpmax)//"rho_model"
      call getmod(vpfile, vsfile, rhofile, vp, vs, rho, fai1, fai2)
!
! Check stability 
!
      call stabil(vp, vs)
!
! Initiate PML related variables
!
      call define_pml(nz, nx,vp,freq0, h,h, dt, npml,&
        a_z, b_z, K_z, a_x, b_x, K_x,&
        a_z_half, b_z_half, K_z_half, a_x_half, b_x_half, K_x_half)
! Compute response functions
      call getfai_vsc(vp,vs,rho,qp,qs1,fai1,fai2,vpe,vse)
!
! Define source function
!
      call series(seri, 1, nt, dt, nt, freq0) ! Ricker

!   
! Output source wavelet
!
      if(myid .eq. 0) then
         filename = datapath(idpmin:idpmax)//"source.bin"
         call centra(filename,imin,imax)
         open(12, file=filename(imin:imax), access='direct',&
             form='unformatted', recl=4*nt )
         write(12, rec=1) (seri(i), i=1, nt)
         close(12)
      endif
!
!
      !loop over shot number
      do isr = myid+1 , nsr, numprcs
!
         call ficnum(isr, cisr)
!
         !define source location in the grid.
         nxsr=xsrg(isr)
         nzsr=zsrg(isr)
!
         !output some info to screen
         write(0,*) "=============="
         write(0,*) "shot number: ", isr
         write(0,*) "at x, z", nxsr, nzsr    
!
         !open f110* for outputting horizontal component
         filename = datapath(idpmin:idpmax)//"f110"//cisr
         call centra(filename, imin, imax)
         open(ifoh, file=filename(imin:imax), access='direct',&
             form='unformatted', recl=mac_len*nx)
!
         !open f120* for outputting vertical component
         filename = datapath(idpmin:idpmax)//"f120"//cisr
         call centra(filename, imin, imax)
         open(ifov, file=filename(imin:imax), access='direct',&
             form='unformatted', recl=mac_len*nx)
!
! Zero recursive arrays
        qs2=real0
        qs1=real0
        qp=real0
        divergv=real0
        sxx=real0
        szz=real0
        szx=real0
        spp=real0
        vx=real0
        vz=real0
        vpx=real0
        vpz=real0
        vsx=real0
        vsz=real0
        memory_dvx_dx=real0
        memory_dvz_dz=real0
        memory_dvx_dz=real0
        memory_dvz_dx=real0
        memory_dsigmazz_dz=real0
        memory_dsigmaxx_dx=real0
        memory_dsigmazx_dx=real0
        memory_dsigmazx_dz=real0
!
! Write snapshots in files
        if (snap .eq. 1) then
                 filename = datapath(idpmin:idpmax)//"s010"//cisr
                 call centra(filename, imin, imax)
                 open(31, file=filename(imin:imax), access='direct',&
                     form='unformatted', recl=mac_len*nx*nz)
                 filename = datapath(idpmin:idpmax)//"s020"//cisr
                 call centra(filename, imin, imax)
                 open(32, file=filename(imin:imax), access='direct',&
                     form='unformatted', recl=mac_len*nx*nz)
                 filename = datapath(idpmin:idpmax)//"s110"//cisr
                 call centra(filename, imin, imax)
                 open(144, file=filename(imin:imax), access='direct',&
                     form='unformatted', recl=mac_len*nx*nz)
                 filename = datapath(idpmin:idpmax)//"s120"//cisr
                 call centra(filename, imin, imax)
                 open(145, file=filename(imin:imax), access='direct',&
                     form='unformatted', recl=mac_len*nx*nz)
                 filename = datapath(idpmin:idpmax)//"s210"//cisr
                 call centra(filename, imin, imax)
                 open(146, file=filename(imin:imax), access='direct',&
                     form='unformatted', recl=mac_len*nx*nz)
                 filename = datapath(idpmin:idpmax)//"s220"//cisr
                 call centra(filename, imin, imax)
                 open(147, file=filename(imin:imax), access='direct',&
                     form='unformatted', recl=mac_len*nx*nz)
        endif

         isnap=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !forward in time to extrapolate the wavefield
         do it = 1, nt
!
            if(mod(it, itsnap) .eq. 0 ) then
                 if (snap .eq. 1) then
                      write(31,rec=isnap) ((vx(i,j), i=1,nz),j=1,nx)
                      write(32,rec=isnap) ((vz(i,j), i=1,nz),j=1,nx)
                      write(144,rec=isnap) ((vpx(i,j), i=1,nz),j=1,nx)
                      write(145,rec=isnap) ((vpz(i,j), i=1,nz),j=1,nx)
                      write(146,rec=isnap) ((vsx(i,j), i=1,nz),j=1,nx)
                      write(147,rec=isnap) ((vsz(i,j), i=1,nz),j=1,nx)

             endif
                isnap=isnap+1
            endif
            if(mod(it, itshow) .eq. 0 ) then
               write(0,*) "myid, isr, it = ", myid, isr, it
            endif
!
          
            !add source
            amp = seri(it)
               ! explosive source
            if (isource .eq. 1) then
                    sxx(nzsr,nxsr)=sxx(nzsr,nxsr)-amp
                    szz(nzsr,nxsr)=szz(nzsr,nxsr)-amp
            else
                ! P-S mixed source
                call s_src(vx,vz,nzsr, nxsr,amp/4)
                call p_src(vx,vz,nzsr, nxsr,amp)
            endif

!
            !**********************!
            !extrapolate wavefield !
            !**********************!
!
            !solve equations of motions
            call eqmot_vis_ps(vp,vs,rho,vx,vz,sxx,szz,szx,&
                       vpx,vpz,vsx,vsz,spp,&
                       memory_dsigmaxx_dx,memory_dsigmazz_dz,&
                       memory_dsigmazx_dx,memory_dsigmazx_dz,&
                       a_z,b_z,k_z,a_x,b_x,k_x,&
                       a_z_half, b_z_half,K_z_half,a_x_half,b_x_half, &
                       K_x_half)



!
! compute memory variables
!
            call memory_ps(vx,vz,memory_dvx_dx,memory_dvz_dz,&
                        memory_dvx_dz,memory_dvz_dx,k_x,k_z,k_x_half,&
                        k_z_half,memory_qs1,memory_qs2,fai1,fai2,qp,&
                        qs1,qs2,divergv,memory_xz,memory_qp)


!
!solve generalized Hooke's relationship
!

            call hooke_vis_ps(vp,vs,rho,vx,vz,sxx,szz,szx,spp,&
                       qp,qs1,qs2,divergv,&
                       memory_dvx_dx,memory_dvz_dz,&
                       memory_dvz_dx,memory_dvx_dz,&
                       a_z,b_z,k_z,a_x,b_x,k_x,&
                       a_z_half, b_z_half,K_z_half,a_x_half,b_x_half,&
                       K_x_half)

!

            !output and save the synthetic observed data to file "f110#".
            rec_z=nzsr
            write(ifoh, rec=it)(vx(rec_z,ix),ix=1,nx)
            write(ifov, rec=it)(vz(rec_z,ix),ix=1,nx)

!
         !end loop over time
         enddo
!
         !close the i/o
         close(ifoh)
         close(ifov)
         close(31)
         close(32)
         close(144)
         close(145)
         close(146)
         close(147)
!
      !end shot loop
      enddo
!
      !synchronize again and get the elapse time
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endtime = MPI_WTIME()  
!
      !write some info. to screen
      if(myid .eq. 0) then
         write(0,*) "time = ", (endtime-starttime)/60.0,"minutes"
         write(0,*) "That's all for lvob2d!"
      endif
!
 301  format('process', i4, ' on ', a)
      call MPI_FINALIZE(ierr)
!
      end
!
!=============================================================================
!

