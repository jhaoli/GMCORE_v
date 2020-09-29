module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use operators_mod
  use reduce_mod

  implicit none

  private

  public damp_init
  public damp_final
  public latlon_damp_lon_limiter
  public latlon_damp_lat_limiter
  public latlon_damp_cell_limiter

  real(r8), allocatable :: cx_full_lat(:,:,:)
  real(r8), allocatable :: cy_full_lat(:,:,:)
  real(r8), allocatable :: cx_half_lat(:,:,:)
  real(r8), allocatable :: cy_half_lat(:,:,:)
contains

  subroutine damp_init()

    integer j, jr, k, r

    call damp_final()

    allocate(cx_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cx_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cx_full_lat(j,k,4) = (global_mesh%full_cos_lat(j) * global_mesh%dlon * 0.5_r8)**4 / dt_in_seconds
          cy_full_lat(j,k,4) = 1.0 / ((1.0_r8 + global_mesh%full_cos_lat(j)**2) / global_mesh%full_cos_lat(j)**2 * 4.0_r8 / &
                          global_mesh%dlat(j)**2 + 16.0_r8 / global_mesh%dlat(j)**4) / dt_in_seconds * 0.01
        end do
      end do
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cx_half_lat(j,k,4) = (global_mesh%half_cos_lat(j) * global_mesh%dlon * 0.5_r8)**4 / dt_in_seconds
          cy_half_lat(j,k,4) = 1.0 / ((1.0_r8 + global_mesh%half_cos_lat(j)**2) / global_mesh%half_cos_lat(j)**2 * 4.0_r8 / &
                          global_mesh%dlat(j)**2 + 16.0_r8 / global_mesh%dlat(j)**4) / dt_in_seconds * 0.01
        end do
      end do

  end subroutine damp_init

  subroutine damp_final()

    if (allocated(cx_full_lat)) deallocate(cx_full_lat)
    if (allocated(cy_full_lat)) deallocate(cy_full_lat)
    if (allocated(cx_half_lat)) deallocate(cx_half_lat)
    if (allocated(cy_half_lat)) deallocate(cy_half_lat) 

  end subroutine damp_final

  subroutine latlon_damp_lon_limiter(block, dt, f)
    ! Fourth diffusion with Xue(2000)'s method
    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy, dfdx, dfdy
    integer i, j, k, ip, cyc

    mesh => block%mesh

    gx   => block%latlon_damp_lon_gx
    gy   => block%latlon_damp_lon_gy
    dfdx => block%latlon_damp_lon_dfdx
    dfdy => block%latlon_damp_lon_dfdy
    cycle_loop: do cyc = 1, polar_damp_cycles
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        ! Calculate damping flux and the gradient of forecast variable at interfaces.
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gx  (i,j,k) = (-f(i-2,j,k) + 3.0 * f(i-1,j,k) - 3.0 * f(i,j,k) + f(i+1,j,k)) / (mesh%full_cos_lat(j) * mesh%dlon)**3
            dfdx(i,j,k) =  f(i,j,k) - f(i-1,j,k)
          end do
        end do

        if (mesh%has_south_pole()) then
          j = mesh%full_lat_ibeg_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            ip = i + mesh%num_half_lon / 2
            if (ip > mesh%num_half_lon) then
              ip = ip - mesh%num_half_lon 
            end if
#ifdef V_POLE
           
#else 
            f(i,j-2,k) = -f(ip,j,k)
#endif
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            ip = i + mesh%num_half_lon / 2
            if (ip > mesh%num_half_lon) then
              ip = ip - mesh%num_half_lon 
            end if
#ifdef V_POLE
            
#else
            f(i,j+2,k) = -f(ip,j,k)
#endif
          end do
        end if 

        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            gy  (i,j,k) = -(f(i,j,k) - f(i,j-1,k)) / mesh%dlat(j) / mesh%half_cos_lat(j)**2 -          &
                           mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                               &
                           (f(i,j+1,k) - f(i,j,k) - f(i,j-1,k) + f(i,j-2,k)) / (2 * mesh%dlat(j)**2) + &
                           (-f(i,j-2,k) + 3 * f(i,j-1,k) - 3 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**3 
            dfdy(i,j,k) = f(i,j,k) - f(i,j-1,k)
#else
            gy  (i,j,k) = -(f(i,j+1,k) - f(i,j,k)) / mesh%dlat(j) / mesh%half_cos_lat(j)**2 -          &
                           mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                               &
                           (f(i,j+2,k) - f(i,j+1,k) - f(i,j,k) + f(i,j-1,k)) / (2 * mesh%dlat(j)**2) + &
                          (-f(i,j-1,k) + 3 * f(i,j,k) - 3 * f(i,j+1,k) + f(i,j+2,k)) / mesh%dlat(j)**3
            dfdy(i,j,k) = f(i,j+1,k) - f(i,j,k)
#endif
          end do
        end do
        ! Limit damping flux to avoid upgradient (Xue 2000).
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * dfdx(i,j,k)))
          end do
        end do
        call fill_halo(block, gx, full_lon=.true., full_lat=.true., full_lev=.true.)
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * dfdy(i,j,k)))
          end do
        end do
        call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true.)
        ! Damp physical variable at last.
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (                                      &
              (gx(i+1,j,k) - gx(i,j,k)) / mesh%full_cos_lat(j) / mesh%dlon * cx_full_lat(j,k,4) + &
              (gy(i,j+1,k) - gy(i,j,k)) / mesh%dlat(j) * cy_full_lat(j,k,4)                       &
              )
#else
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (                                      &
              (gx(i+1,j,k) - gx(i,j,k)) / mesh%full_cos_lat(j) / mesh%dlon * cx_full_lat(j,k,4) + &
              (gy(i,j,k) - gy(i,j-1,k)) / mesh%dlat(j) * cy_full_lat(j,k,4)                       &
            ) 
#endif
          end do
        end do
      end do
      call fill_halo(block, f, full_lon=.false., full_lat=.true., full_lev=.true.)
    end do cycle_loop
    
  end subroutine latlon_damp_lon_limiter

  subroutine latlon_damp_lat_limiter(block, dt, f)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy, dfdx, dfdy
    integer i, j, k, cyc, ip

    mesh => block%mesh

    gx   => block%latlon_damp_lat_gx
    gy   => block%latlon_damp_lat_gy
    dfdx => block%latlon_damp_lat_dfdx
    dfdy => block%latlon_damp_lat_dfdy
    cycle_loop: do cyc = 1, polar_damp_cycles
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        ! Calculate damping flux and the gradient of forecast variable at interfaces.
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gx  (i,j,k) = (-f(i-1,j,k) + 3 * f(i,j,k) - 3 * f(i+1,j,k) + f(i+2,j,k)) / (mesh%half_cos_lat(j) * mesh%dlon)**3
            dfdx(i,j,k) = f(i+1,j,k) - f(i,j,k)
          end do
        end do

        if (mesh%has_south_pole()) then
          j = mesh%half_lat_ibeg_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip = i + mesh%num_full_lon / 2
            if (ip > mesh%num_full_lon) then
              ip = ip - mesh%num_full_lon
            end if
#ifdef V_POLE
              
#else
              f(i,j-1,k) = -f(ip,j  ,k)
#endif
          end do
        end if
        if (mesh%has_north_pole()) then 
          j = mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip = i + mesh%num_full_lon / 2
            if (ip > mesh%num_full_lon) then
              ip = ip - mesh%num_full_lon
            end if
#ifdef V_POLE

#else
            f(i,j+1,k) = -f(ip,j  ,k)
#endif    
          end do           
        end if 

        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            gy  (i,j,k) = -(f(i,j+1,k) - f(i,j,k)) / mesh%full_cos_lat(j)**2 / mesh%dlat(j) -          &
                          mesh%full_sin_lat(j) / mesh%full_cos_lat(j) *                                &
                          (f(i,j+2,k) - f(i,j+1,k) - f(i,j,k) + f(i,j-1,k)) / (2 * mesh%dlat(j)**2) +  &
                          (-f(i,j-1,k) + 3 * f(i,j,k) - 3 * f(i,j+1,k) + f(i,j+2,k)) / mesh%dlat(j)**3 
            dfdy(i,j,k) = f(i,j+1,k) - f(i,j,k)
#else
            gy  (i,j,k) = -(f(i,j,k) - f(i,j-1,k)) / mesh%full_cos_lat(j)**2 / mesh%dlat(j) -          &
                           mesh%full_sin_lat(j) / mesh%full_cos_lat(j) *                               &
                           (f(i,j+1,k) - f(i,j,k) - f(i,j-1,k) + f(i,j-2,k)) / (2 * mesh%dlat(j)**2) + &
                          (-f(i,j-2,k) + 3 * f(i,j-1,k) - 3 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**3
            dfdy(i,j,k) = f(i,j,k) - f(i,j-1,k)
#endif
          end do
        end do
        ! Limit damping flux to avoid upgradient (Xue 2000).
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * dfdx(i,j,k)))
          end do
        end do
        call fill_halo(block, gx, full_lon=.false., full_lat=.false., full_lev=.true.)
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * dfdy(i,j,k)))
          end do
        end do
        call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true.)
        ! Damp physical variable at last.
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (                                      &
              (gx(i,j,k) - gx(i-1,j,k)) / mesh%half_cos_lat(j) / mesh%dlon * cx_half_lat(j,k,4) + &
              (gy(i,j,k) - gy(i,j-1,k)) / mesh%dlat(j) * cy_half_lat(j,k,4)                       &
              )
#else
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (                                      &
              (gx(i,j,k) - gx(i-1,j,k)) / mesh%half_cos_lat(j) / mesh%dlon * cx_half_lat(j,k,4) + &
              (gy(i,j+1,k) - gy(i,j,k)) / mesh%dlat(j) * cy_half_lat(j,k,4)                       &
            ) 
#endif
          end do
        end do
      end do
      call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)
    end do cycle_loop

  end subroutine latlon_damp_lat_limiter

  subroutine latlon_damp_cell_limiter(block, dt, f)
    ! Fourth diffusion with Xue(2000)'s method
    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy, dfdx, dfdy
    integer i, j, k, ip, cyc

    mesh => block%mesh

    gx   => block%latlon_damp_cell_gx
    gy   => block%latlon_damp_cell_gy
    dfdx => block%latlon_damp_cell_dfdx
    dfdy => block%latlon_damp_cell_dfdy
    cycle_loop: do cyc = 1, polar_damp_cycles
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        ! Calculate damping flux and the gradient of forecast variable at interfaces.
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gx  (i,j,k) = (-f(i-1,j,k) + 3.0 * f(i,j,k) - 3.0 * f(i+1,j,k) + f(i+2,j,k)) / (mesh%full_cos_lat(j) * mesh%dlon)**3
            dfdx(i,j,k) =  f(i+1,j,k) - f(i,j,k)
          end do
        end do

        if (mesh%has_south_pole()) then
          j = mesh%full_lat_ibeg_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip = i + mesh%num_full_lon / 2
            if (ip > mesh%num_full_lon) then
              ip = ip - mesh%num_full_lon 
            end if
#ifdef V_POLE
             
#else 
            f(i,j-2,k) = f(ip,j,k)
#endif
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip = i + mesh%num_full_lon / 2
            if (ip > mesh%num_full_lon) then
              ip = ip - mesh%num_full_lon 
            end if
#ifdef V_POLE
            
#else
            f(i,j+2,k) = f(ip,j,k)
#endif
          end do
        end if 

        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            gy  (i,j,k) = -(f(i,j,k) - f(i,j-1,k)) / mesh%dlat(j) / mesh%half_cos_lat(j)**2 -          &
                           mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                               &
                           (f(i,j+1,k) - f(i,j,k) - f(i,j-1,k) + f(i,j-2,k)) / (2 * mesh%dlat(j)**2) + &
                           (-f(i,j-2,k) + 3 * f(i,j-1,k) - 3 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**3 
            dfdy(i,j,k) = f(i,j,k) - f(i,j-1,k)
#else
            gy  (i,j,k) = -(f(i,j+1,k) - f(i,j,k)) / mesh%dlat(j) / mesh%half_cos_lat(j)**2 -          &
                           mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                               &
                           (f(i,j+2,k) - f(i,j+1,k) - f(i,j,k) + f(i,j-1,k)) / (2 * mesh%dlat(j)**2) + &
                          (-f(i,j-1,k) + 3 * f(i,j,k) - 3 * f(i,j+1,k) + f(i,j+2,k)) / mesh%dlat(j)**3
            dfdy(i,j,k) = f(i,j+1,k) - f(i,j,k)
#endif
          end do
        end do
        ! Limit damping flux to avoid upgradient (Xue 2000).
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * dfdx(i,j,k)))
          end do
        end do
        call fill_halo(block, gx, full_lon=.true., full_lat=.true., full_lev=.true.)
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * dfdy(i,j,k)))
          end do
        end do
        call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true.)
        ! Damp physical variable at last.
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            ! f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (                                      &
            !   (gx(i,j,k) - gx(i-1,j,k)) / mesh%full_cos_lat(j) / mesh%dlon * cx_full_lat(j,k,4) + &
            !   (gy(i,j,k) - gy(i,j-1,k)) / mesh%dlat(j) * cy_full_lat(j,k,4)                       &
            !   )
#else
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (                                      &
              (gx(i,j,k) - gx(i-1,j,k)) / mesh%full_cos_lat(j) / mesh%dlon * cx_full_lat(j,k,4) + &
              (gy(i,j,k) - gy(i,j-1,k)) / mesh%dlat(j) * cy_full_lat(j,k,4)                       &
            ) 
#endif
          end do
        end do
      end do
      call fill_halo(block, f, full_lon=.true., full_lat=.true., full_lev=.true.)
    end do cycle_loop

  end subroutine latlon_damp_cell_limiter

end module damp_mod
