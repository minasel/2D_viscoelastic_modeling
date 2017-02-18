      integer nx,nz
      integer isnap
      real, allocatable :: h(:,:,:),v(:,:,:)
      character*128 filename
      character*128 inputname1
      character*128 inputname2
      character*4 isnap_char 

      read (*,*) inputname1
      read (*,*) inputname2
      read (*,*) nz, nx, nsnap
      allocate (h(nz,nx,nsnap))
      allocate (v(nz,nx,nsnap))
ccccccccccccccc Snap cccccccccccccccccccccc
      open(8,file=inputname1,access='direct',recl=4*nx*nz)
      open(9,file=inputname2,access='direct',recl=4*nx*nz)
      do isnap=1,nsnap
        read(8,rec=isnap) ((h(i,j,isnap),i=1,nz), j=1,nx)
        read(9,rec=isnap) ((v(i,j,isnap),i=1,nz), j=1,nx)
      enddo
        !write(11,rec=1) nx

      do isnap=1,nsnap
        call ficnum(isnap,isnap_char)
        filename="snap_v_"//isnap_char
        open(isnap+10,file=filename,access='direct',recl=4*nx*nz)
        write(isnap+10,rec=1) ((v(i,j,isnap),i=1,nz), j=1,nx)
        close (isnap+10)
      enddo

      do isnap=1,nsnap
        call ficnum(isnap,isnap_char)
        filename="snap_h_"//isnap_char
        open(isnap,file=filename,access='direct',recl=4*nx*nz)
        write(isnap,rec=1) ((h(i,j,isnap),i=1,nz), j=1,nx)
        close(isnap)
      enddo

        close(8)
        close(9)
        return
        end
