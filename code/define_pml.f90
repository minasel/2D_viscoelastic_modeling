! subroutine to define CPML boundaries
! From SEISMIC-CPML package (See introduction of documentation)
! Copyright see license CeCILL.text in ../doc/ directory

 SUBROUTINE define_pml(nz, nx, vp, f0, dz, dx, dt,NPOINTS_PML,&
        a_z,b_z, K_z, a_x, b_x, K_x,&
        a_z_half, b_z_half, K_z_half, a_x_half, b_x_half, K_x_half)

        implicit none

        ! input
        integer :: nz, nx

        real :: dz, dx, dt, f0
        real :: vp(-4:nz+4,-4:nx+4)

        ! output
        real, dimension(nz) :: a_z, b_z, K_z
        real, dimension(nx) :: a_x, b_x, K_x
        real, dimension(nz) :: a_z_half, b_z_half, K_z_half
        real, dimension(nx) :: a_x_half, b_x_half, K_x_half

        ! flags to add PML layers to the edges of the grid
        logical, parameter :: USE_PML_ZMIN = .true.
        logical, parameter :: USE_PML_ZMAX = .true.
        logical, parameter :: USE_PML_XMIN = .true.
        logical, parameter :: USE_PML_XMAX = .true.

        ! thickness of the PML layer in grid points
        !integer, parameter :: NPOINTS_PML = 15
        integer :: NPOINTS_PML

        ! value of PI
        real, parameter :: PI = 4. / atan(1.)

        ! zero
        real, parameter :: real0 = 0.d0

        ! power to compute d0 profile
        real, parameter :: NPOWER = 2.d0

        real, parameter :: K_MAX_PML = 1.d0 ! from Gedney page 8.11
        !from Festa and Vilotte:
!!!        real, parameter :: ALPHA_MAX_PML = 2.e0*PI*(f0/2.e0)
        real  ALPHA_MAX_PML

        ! 1D arrays for the damping profiles
        real, dimension(nz) :: d_z, alpha_z, d_z_half, alpha_z_half
        real, dimension(nx) :: d_x, alpha_x, d_x_half, alpha_x_half

        real :: thickness_PML_z, thickness_PML_x
        real :: zoriginleft, zoriginright, xoriginbottom, xorigintop
        real :: Rcoef, d0_z(nz,2), d0_x(nx,2), zval, xval
        real :: abscissa_in_PML, abscissa_normalized

        integer :: iz, ix

!--- define profile of absorption in PML region
        ALPHA_MAX_PML = 2.e0*PI*(f0/2.e0)
! thickness of the PML layer in meters
  thickness_PML_z = NPOINTS_PML * dz
  thickness_PML_x = NPOINTS_PML * dx

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 1e-5

! check that NPOWER is okay
  if(NPOWER < 1) stop 'NPOWER must be greater than 1'

!       compute d0 from INRIA report section 6.1 
!       http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf

  do iz = 1,nz
     !left:
     d0_z(iz,1) = -(NPOWER+1) * vp(iz,1) * log(Rcoef) / (2.e0 * thickness_PML_z)

     !right:
     d0_z(iz,2) = -(NPOWER+1) * vp(iz,nx) * log(Rcoef) / (2.e0 * thickness_PML_z)
  enddo

  do ix = 1,nx
     !bottom:
     d0_x(ix,1) = -(NPOWER+1) * vp(nz,ix) * log(Rcoef) / (2.e0 * thickness_PML_x)

     !top:
     d0_x(ix,2) = -(NPOWER+1) * vp(1,ix) * log(Rcoef) / (2.e0 * thickness_PML_x)
  enddo

!!!  print *,'d0_z = ',d0_z
!!!  print *,'d0_x = ',d0_x
!!!  print *

  d_z(:) = real0
  d_z_half(:) = real0
  K_z(:) = 1.e0
  K_z_half(:) = 1.e0
  alpha_z(:) = real0
  alpha_z_half(:) = real0
  a_z(:) = real0
  a_z_half(:) = real0

  d_x(:) = real0
  d_x_half(:) = real0
  K_x(:) = 1.e0
  K_x_half(:) = 1.e0
  alpha_x(:) = real0
  alpha_x_half(:) = real0
  a_x(:) = real0
  a_x_half(:) = real0
! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  zoriginleft = thickness_PML_z
  zoriginright = (nz-1)*dz - thickness_PML_z

  do iz = 1,nz

! abscissa of current grid point along the damping profile
    zval = dz * real(iz-1)

!---------- left edge
    if(USE_PML_ZMIN) then

! define damping profile at the grid points
      abscissa_in_PML = zoriginleft - zval
      if(abscissa_in_PML >= real0) then
!        print*,iz
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z(iz) = d0_z(iz,1) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_z(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
      endif
! define damping profile at half the grid points
      abscissa_in_PML = zoriginleft - (zval + dz/2.e0)
      if(abscissa_in_PML >= real0) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z_half(iz) = d0_z(iz,1) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z_half(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_z_half(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0* ALPHA_MAX_PML
      endif


    endif

!---------- right edge
    if(USE_PML_ZMAX) then

! define damping profile at the grid points
      abscissa_in_PML = zval - zoriginright
      if(abscissa_in_PML >= real0) then
!        print*,iz
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z(iz) = d0_z(iz,2) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_z(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = zval + dz/2.e0 - zoriginright
      if(abscissa_in_PML >= real0) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z_half(iz) = d0_z(iz,2) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z_half(iz) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_z_half(iz) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0* ALPHA_MAX_PML
      endif

    endif

! just in case, for -5 at the end
    if(alpha_z(iz) < real0) alpha_z(iz) = real0
    if(alpha_z_half(iz) < real0) alpha_z_half(iz) = real0

    b_z(iz) = exp(- (d_z(iz) / K_z(iz) + alpha_z(iz)) * dt)
    b_z_half(iz) = exp(- (d_z_half(iz) / K_z_half(iz) + alpha_z_half(iz)) * dt)

! this to avoid division by zero outside the PML
    if(abs(d_z(iz)) > 1.e-6) a_z(iz) = d_z(iz) * (b_z(iz) - 1.e0) / (K_z(iz) * (d_z(iz) + K_z(iz) * alpha_z(iz)))
    if(abs(d_z_half(iz)) > 1.e-6) a_z_half(iz) = d_z_half(iz) * &
      (b_z_half(iz) - 1.e0) / (K_z_half(iz) * (d_z_half(iz) + K_z_half(iz) * alpha_z_half(iz)))

  enddo

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginbottom = thickness_PML_x
  xorigintop = (nx-1)*dx - thickness_PML_x

  do ix = 1,nx

! abscissa of current grid point along the damping profile
    xval = dx * real(ix-1)

!---------- bottom edge
    if(USE_PML_XMIN) then

! define damping profile at the grid points
      abscissa_in_PML = xoriginbottom - xval
      if(abscissa_in_PML >= real0) then
!        print*,ix
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(ix) = d0_x(ix,1) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_x(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
      endif
! define damping profile at half the grid points
      abscissa_in_PML = xoriginbottom - (xval + dx/2.e0)
      if(abscissa_in_PML >= real0) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(ix) = d0_x(ix,1) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x_half(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_x_half(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0* ALPHA_MAX_PML
      endif


    endif

!---------- top edge
    if(USE_PML_XMAX) then

! define damping profile at the grid points
      abscissa_in_PML = xval - xorigintop
      if(abscissa_in_PML >= real0) then
!        print*,ix
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(ix) = d0_x(ix,2) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_x(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
      endif
! define damping profile at half the grid points
      abscissa_in_PML = xval + dx/2.e0 - xorigintop
      if(abscissa_in_PML >= real0) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(ix) = d0_x(ix,2) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x_half(ix) = 1.e0 + (K_MAX_PML - 1.e0) * abscissa_normalized**NPOWER
        alpha_x_half(ix) = ALPHA_MAX_PML * (1.e0 - abscissa_normalized) + 0.1e0 * ALPHA_MAX_PML
      endif

    endif

    b_x(ix) = exp(- (d_x(ix) / K_x(ix) + alpha_x(ix)) * dt)
    b_x_half(ix) = exp(- (d_x_half(ix) / K_x_half(ix) + alpha_x_half(ix)) * dt)

! this to avoid division by zero outside the PML
    if(abs(d_x(ix)) > 1.e-6) a_x(ix) = d_x(ix) * (b_x(ix) - 1.e0) / (K_x(ix) * (d_x(ix) + K_x(ix) * alpha_x(ix)))
    if(abs(d_x_half(ix)) > 1.e-6) a_x_half(ix) = d_x_half(ix) * &
      (b_x_half(ix) - 1.e0) / (K_x_half(ix) * (d_x_half(ix) + K_x_half(ix) *alpha_x_half(ix)))


  enddo
        RETURN
        END SUBROUTINE define_pml
!***********************************************************************
!***********************************************************************


