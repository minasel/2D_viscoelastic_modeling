!Model
	real, allocatable::&
		zsrg(:),xsrg(:),seri(:),&
		vp(:,:), vs(:,:),rho(:,:),vpe(:,:),vse(:,:),&
		qp(:,:), qs1(:,:),qs2(:,:),divergv(:,:),&
		sxx(:,:), szz(:,:),szx(:,:),spp(:,:),sss(:,:),&
		sxx_q(:,:), szz_q(:,:), szx_q(:,:),&
		vx(:,:), vz(:,:),vpx(:,:),vpz(:,:),&
		vsx(:,:), vsz(:,:),pwave(:,:),swave(:,:),&
		fai1(:,:,:), fai2(:,:,:),&
		memory_qp(:,:,:),memory_qs1(:,:,:),memory_qs2(:,:,:),memory_xz(:,:,:),&
		vpx1(:,:),vpz1(:,:),vsx1(:,:), vsz1(:,:)
	real, allocatable::&
		kx_2(:,:),kz_2(:,:),kx(:),kz(:)

!PML Parameters
	real, allocatable::&
      		memory_dvx_dx(:,:), memory_dvz_dz(:,:),&
      		memory_dvx_dz(:,:), memory_dvz_dx(:,:),&
      		memory_dsigmaxx_dx(:,:), memory_dsigmazz_dz(:,:),&
      		memory_dsigmazx_dx(:,:), memory_dsigmazx_dz(:,:),&
      		a_z(:),b_z(:),k_z(:), a_x(:),b_x(:),k_x(:),&
      		a_z_half(:),b_z_half(:),k_z_half(:),&
	 	a_x_half(:),b_x_half(:),k_x_half(:)
! For visco
	real, allocatable::&
      		memory_dvx_dx_q(:,:), memory_dvz_dz_q(:,:),&
      		memory_dvx_dz_q(:,:), memory_dvz_dx_q(:,:),&
      		memory_dsigmaxx_dx_q(:,:), memory_dsigmazz_dz_q(:,:),&
      		memory_dsigmazx_dx_q(:,:), memory_dsigmazx_dz_q(:,:)



