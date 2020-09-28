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
  public polar_damp
  public div_damp
  public latlon_damp_lon_limiter
  public latlon_damp_lat_limiter
  public latlon_damp_cell_limiter

  integer, parameter :: diff_halo_width(2:8) = [1, 2, 2, 3, 3, 4, 4]

  real(r8), target :: diff_weights(9,2:8) = reshape([  &
     1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
    -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
     1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
    -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
     1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
    -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
     1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
  ], [9, 7])

  logical , allocatable :: zonal_damp_on_full_lat(:)
  logical , allocatable :: zonal_damp_on_half_lat(:)
  real(r8), allocatable :: cd_full_lat(:,:)
  real(r8), allocatable :: cd_half_lat(:,:)

  real(r8), allocatable :: cx_full_lat(:,:,:)
  real(r8), allocatable :: cy_full_lat(:,:,:)
  real(r8), allocatable :: cx_half_lat(:,:,:)
  real(r8), allocatable :: cy_half_lat(:,:,:)
contains

  subroutine damp_init()

    integer j, jr, k, r
    real(r8) c0

    call damp_final()

    allocate(zonal_damp_on_full_lat(global_mesh%num_full_lat)); zonal_damp_on_full_lat = .false.
    allocate(zonal_damp_on_half_lat(global_mesh%num_half_lat)); zonal_damp_on_half_lat = .false.

    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%full_lat_iend_no_pole - j + 1
      end if
      if (reduce_factors(jr) > 1) then
        zonal_damp_on_full_lat(j) = .true.
      end if
    end do
    do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
      if (global_mesh%half_lat(j) <= 0) then
        jr = j - global_mesh%half_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%half_lat_iend_no_pole - j + 1
      end if
#ifdef V_POLE
      if (reduce_factors(jr) > 1 .or. reduce_factors(jr-1) > 1) then
#else
      if (reduce_factors(jr) > 1 .or. reduce_factors(jr+1) > 1) then
#endif
        zonal_damp_on_half_lat(j) = .true.
      end if
    end do

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

    allocate(cd_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cd_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cd_full_lat(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cd_half_lat(j,k) = div_damp_coef2 * global_mesh%half_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case (4)
      r = 2

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cd_full_lat(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cd_half_lat(j,k) = div_damp_coef4 * global_mesh%half_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_string(div_damp_order)) // '!')
    end select

  end subroutine damp_init

  subroutine damp_final()

    if (allocated(cx_full_lat)) deallocate(cx_full_lat)
    if (allocated(cy_full_lat)) deallocate(cy_full_lat)
    if (allocated(cx_half_lat)) deallocate(cx_half_lat)
    if (allocated(cy_half_lat)) deallocate(cy_half_lat) 
    if (allocated(cd_full_lat)) deallocate(cd_full_lat)
    if (allocated(cd_half_lat)) deallocate(cd_half_lat)

  end subroutine damp_final
  
  subroutine polar_damp(block, dt, state)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    call zonal_damp(block, polar_damp_order, dt, state%u, &
                    block%mesh%half_lon_lb, block%mesh%half_lon_ub, &
                    block%mesh%full_lat_lb, block%mesh%full_lat_ub, &
                    block%mesh%full_lev_lb, block%mesh%full_lev_ub, &
                    zonal_damp_on_full_lat)
    call zonal_damp(block, polar_damp_order, dt, state%v, &
                    block%mesh%full_lon_lb, block%mesh%full_lon_ub, &
                    block%mesh%half_lat_lb, block%mesh%half_lat_ub, &
                    block%mesh%full_lev_lb, block%mesh%full_lev_ub, &
                    zonal_damp_on_half_lat)
    if (baroclinic) then
      call zonal_damp(block, polar_damp_order, dt, state%pt, &
                      block%mesh%full_lon_lb, block%mesh%full_lon_ub, &
                      block%mesh%full_lat_lb, block%mesh%full_lat_ub, &
                      block%mesh%full_lev_lb, block%mesh%full_lev_ub, &
                      zonal_damp_on_full_lat)
    end if

  end subroutine polar_damp

  subroutine zonal_damp(block, order, dt, f, is, ie, js, je, ks, ke, on)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(is:ie,js:je,ks:ke)
    integer, intent(in) :: is, ie, js, je, ks, ke
    logical, intent(in) :: on(:)

    type(mesh_type), pointer :: mesh
    real(r8) c, gx(is+block%mesh%lon_halo_width:ie-block%mesh%lon_halo_width+1)
    integer i, j, k, o, ns
    real(r8), pointer :: w(:)

    o = order
    c = 0.5_r8**o / dt
    if (o == 2) then
      ns = diff_halo_width(2)
      w => diff_weights(:,2)
    else
      ns = diff_halo_width(o-1)
      w => diff_weights(:,o-1)
    end if

    mesh => block%mesh
    do k = ks, ke
      do j = js + mesh%lat_halo_width, je - mesh%lat_halo_width
        if (on(j)) then
          ! Calculate damping flux at interfaces.
          do i = lbound(gx, 1), ubound(gx, 1)
            gx(i) = sum(f(i-ns:i+ns-1,j,k) * w(:2*ns))
          end do
          if (o > 2) then
            ! Limit damping flux to avoid upgradient (Xue 2000).
            gx = gx * (-1)**(o / 2)
            do i = lbound(gx, 1), ubound(gx, 1)
              gx(i) = gx(i) * max(0.0_r8, sign(1.0_r8, -gx(i) * (f(i,j,k) - f(i-1,j,k))))
            end do
          end if
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            f(i,j,k) = f(i,j,k) - c * dt * (gx(i+1) - gx(i))
          end do
          call fill_zonal_halo(block, mesh%lon_halo_width, f(:,j,k))
        end if
      end do
    end do

  end subroutine zonal_damp

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
  
  subroutine div_damp(block, dt, old_state, new_state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    type(mesh_type), pointer :: mesh
    integer i, j, k, cyc

    if (baroclinic) then
      mesh => old_state%mesh

      select case (div_damp_order)
      case (2)
        new_state%div = old_state%div
        cycle_loop: do cyc = 1, div_damp_cycles
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                new_state%u(i,j,k) = new_state%u(i,j,k) + dt / div_damp_cycles * cd_full_lat(j,k) * ( &
                  new_state%div(i+1,j,k) - new_state%div(i,j,k)) / mesh%de_lon(j)
              end do
            end do
          end do
          call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                new_state%v(i,j,k) = new_state%v(i,j,k) + dt / div_damp_cycles * cd_half_lat(j,k) * ( &
                  new_state%div(i,j,k) - new_state%div(i,j-1,k)) / mesh%de_lat(j)
#else
                new_state%v(i,j,k) = new_state%v(i,j,k) + dt / div_damp_cycles * cd_half_lat(j,k) * ( &
                  new_state%div(i,j+1,k) - new_state%div(i,j,k)) / mesh%de_lat(j)
#endif
              end do
            end do
          end do
          call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

          call calc_div(block, new_state)
        end do cycle_loop
      case (4)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              new_state%u(i,j,k) = new_state%u(i,j,k) - dt * cd_full_lat(j,k) * ( &
                old_state%div2(i+1,j,k) - old_state%div2(i,j,k)) / mesh%de_lon(j)
            end do
          end do
        end do
        call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              new_state%v(i,j,k) = new_state%v(i,j,k) - dt * cd_half_lat(j,k) * ( &
                old_state%div2(i,j,k) - old_state%div2(i,j-1,k)) / mesh%de_lat(j)
#else
              new_state%v(i,j,k) = new_state%v(i,j,k) - dt * cd_half_lat(j,k) * ( &
                old_state%div2(i,j+1,k) - old_state%div2(i,j,k)) / mesh%de_lat(j)
#endif
            end do
          end do
        end do
#ifndef V_POLE
        if (mesh%has_south_pole()) then
          j = mesh%half_lat_ibeg
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%v(i,j,k) = new_state%v(i,j,k) + dt * 3.5d4 * ( &
                old_state%div(i,j+1,k) - old_state%div(i,j,k)) / mesh%de_lat(j)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%half_lat_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%v(i,j,k) = new_state%v(i,j,k) + dt * 3.5d4 * ( &
                old_state%div(i,j+1,k) - old_state%div(i,j,k)) / mesh%de_lat(j)
            end do
          end do
        end if
#endif
        call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end select
    end if

  end subroutine div_damp

end module damp_mod
