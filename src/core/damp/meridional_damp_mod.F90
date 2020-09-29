module meridional_damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public meridional_damp_init
  public meridional_damp_final
  public meridional_damp_on_lon_edge
  public meridional_damp_on_lat_edge
  public meridional_damp_on_cell
  public meridional_damp_on_vtx

  real(r8), allocatable :: cy_full_lat(:,:,:)
  real(r8), allocatable :: cy_half_lat(:,:,:)

contains

  subroutine meridional_damp_init()
  
    integer j, k
    call meridional_damp_final()

    allocate(cy_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
        cy_full_lat(j,k,2) = (global_mesh%dlat(j) / 2.0_r8)**2 / dt_in_seconds * 0.01
        cy_full_lat(j,k,4) = 1.0 / ((1.0_r8 + global_mesh%full_cos_lat(j)**2) / global_mesh%full_cos_lat(j)**2 * 4.0_r8 / &
                          global_mesh%dlat(j)**2 + 16.0_r8 / global_mesh%dlat(j)**4) / dt_in_seconds * 0.01
      end do
      do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
        cy_half_lat(j,k,2) = (global_mesh%dlat(j) / 2.0_r8)**2 / dt_in_seconds * 0.01
        cy_half_lat(j,k,4) = 1.0 / ((1.0_r8 + global_mesh%half_cos_lat(j)**2) / global_mesh%half_cos_lat(j)**2 * 4.0_r8 / &
                          global_mesh%dlat(j)**2 + 16.0_r8 / global_mesh%dlat(j)**4) / dt_in_seconds * 0.01
      end do
    end do
    
  end subroutine meridional_damp_init

  subroutine meridional_damp_final()
    
    if (allocated(cy_full_lat)) deallocate(cy_full_lat)
    if (allocated(cy_half_lat)) deallocate(cy_half_lat)

  end subroutine meridional_damp_final

  subroutine meridional_damp_on_lon_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gy
    integer i, io, j, k

    mesh => block%mesh
    gy   => block%meridional_damp_lon_gy

    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i + global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j-1,:) = -f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i + global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j+1,:) = -f(io,j-1,:)
      end do
    end if
    
    if(order == 4) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        ! Calculate damping flux 
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gy(i,j,k) = -(f(i,j+1,k) - f(i,j,k)) / mesh%dlat(j) / mesh%half_cos_lat(j)**2 -          &
                         mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                               &
                         (f(i,j+2,k) - f(i,j+1,k) - f(i,j,k) + f(i,j-1,k)) / (2 * mesh%dlat(j)**2) + &
                         (-f(i,j-1,k) + 3 * f(i,j,k) - 3 * f(i,j+1,k) + f(i,j+2,k)) / mesh%dlat(j)**3
          end do
        end do

        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j+1,k) - f(i,j,k))))
          end do
        end do
        call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true.)
        ! Damp physical variables at last.
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            f(i,j,k) = f(i,j,k) - dt * (gy(i,j,k) - gy(i,j-1,k)) / mesh%dlat(j) * cy_full_lat(j,k,4)
          end do
        end do
      end do
    else if (order == 2) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
!            f(i,j,k) = f(i,j,k) + dt * (-mesh%full_sin_lat(j) / mesh%full_cos_lat(j) * (f(i,j+1,k) - f(i,j-1,k)) / &
!                    (2.d0 * mesh%dlat(j)) + (f(i,j-1,k) - 2.d0 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**2) * cy_full_lat(j,k,2)
            f(i,j,k) = f(i,j,k) + dt / mesh%full_cos_lat(j) / mesh%dlat(j)**2 *&
                              (mesh%half_cos_lat(j  ) * (f(i,j+1,k) - f(i,j  ,k)) -&
                               mesh%half_cos_lat(j-1) * (f(i,j  ,k) - f(i,j-1,k))) * cy_full_lat(j,k,2)         
          end do
        end do
      end do
    end if

    call fill_halo(block, f, full_lon=.false., full_lat=.true., full_lev=.true.)

  end subroutine meridional_damp_on_lon_edge

  subroutine meridional_damp_on_lat_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gy
    integer i, io, j, k

    mesh => block%mesh 
    gy   => block%meridional_damp_lat_gy

    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        io = i + global_mesh%num_full_lon / 2
        if (io > global_mesh%num_full_lon) io = io - global_mesh%num_full_lon
        f(i,j-1,:) = -f(io,j  ,:)
        f(i,j-2,:) = -f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        io = i + global_mesh%num_full_lon / 2
        if (io > global_mesh%num_full_lon) io = io - global_mesh%num_full_lon
        f(i,j+1,:) = -f(io,j  ,:)
        f(i,j+2,:) = -f(io,j-1,:)
      end do
    end if
    
    if (order == 4) then
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Calculate damping flux
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gy(i,j,k) = -1 / mesh%full_cos_lat(j)**2 * (f(i,j,k) - f(i,j-1,k)) / mesh%dlat(j)  -    &
                      mesh%full_sin_lat(j) / mesh%full_cos_lat(j) *                               &
                      (f(i,j+1,k) - f(i,j,k) - f(i,j-1,k) + f(i,j-2,k)) / (2 * mesh%dlat(j)**2) + &
                      (-f(i,j-2,k) + 3 * f(i,j-1,k) - 3 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**3
        end do
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j,k) - f(i,j-1,k))))
        end do
      end do
      call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true.)
      ! Damp physical variable at last.
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          f(i,j,k) = f(i,j,k) - dt * (gy(i,j+1,k) - gy(i,j,k)) / mesh%dlat(j) * cy_half_lat(j,k,4)
        end do
      end do
    end do
    else if (order == 2) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
!            f(i,j,k) = f(i,j,k) + dt * (-mesh%half_sin_lat(j) / mesh%half_cos_lat(j) * (f(i,j+1,k) - f(i,j-1,k)) /&
!                    (2.d0 * mesh%dlat(j)) + (f(i,j-1,k) - 2.d0 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**2) * cy_half_lat(j,k,2)
            f(i,j,k) = f(i,j,k) + dt / mesh%half_cos_lat(j) / mesh%dlat(j)**2 * &
                          (mesh%full_cos_lat(j+1) * (f(i,j+1,k) - f(i,j  ,k)) - &
                           mesh%full_cos_lat(j  ) * (f(i,j  ,k) - f(i,j-1,k))) * cy_half_lat(j,k,2)
          end do
        end do
      end do
    end if 
    call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine meridional_damp_on_lat_edge

  subroutine meridional_damp_on_cell(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    
    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gy
    integer i, io, j, k

    mesh => block%mesh
    gy   => block%meridional_damp_cell_gy

    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        io = io + mesh%num_full_lon / 2
        if (io > mesh%num_full_lon) io = io - mesh%num_full_lon
        f(i,j-1,:) = f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        io = i + mesh%num_full_lon / 2
        if (io > mesh%num_full_lon) io = io - mesh%num_full_lon
        f(i,j+1,:) = f(i,j-1,:)
      end do
    end if
    
    if (order == 4) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        ! Calculate damping flux
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gy(i,j,k) = -1/ mesh%half_cos_lat(j)**2 * (f(i,j+1,k) - f(i,j,k)) / mesh%dlat(j) -      &
                        mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                               &
                        (f(i,j+2,k) - f(i,j+1,k) - f(i,j,k) + f(i,j-1,k)) / (2 * mesh%dlat(j)**2) + &
                        (-f(i,j-1,k) + 3 * f(i,j,k) - 3 * f(i,j+1,k) + f(i,j+2,k)) / mesh%dlat(j)**3
          end do
        end do
        ! Limit damping flux to avoid upgradient (Xue 2000).
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j+1,k) - f(i,j,k))))
          end do
        end do
        call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true.)
        ! Damp physical variable at last.
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            f(i,j,k) = f(i,j,k) - dt * (gy(i,j,k) - gy(i,j-1,k)) / mesh%dlat(j) * cy_full_lat(j,k,4)
          end do
        end do
      end do
    else if (order == 2) then
      do k = mesh%full_lev_ibeg, mesh%full_lon_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            f(i,j,k) = f(i,j,k) + dt * (mesh%half_cos_lat(j) * (f(i,j+1,k) - f(i,j,k)) -&
                                        mesh%half_cos_lat(j-1) * (f(i,j,k) - f(i,j-1,k))) /&
                                        mesh%full_cos_lat(j) / mesh%dlat(j)**2 * cy_full_lat(j,k,2) 
          end do
        end do
      end do
    end if
    call fill_halo(block, f, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine meridional_damp_on_cell

  subroutine meridional_damp_on_vtx(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gy
    integer i, io, j, k

    mesh => block%mesh
    gy   => block%meridional_damp_vtx_gy

    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i + global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j-1,:) = f(io,j  ,:)
        f(i,j-2,:) = f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i + global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j+1,:) = f(io,j  ,:)
        f(i,j+2,:) = f(io,j-1,:)
      end do
    end if
    
    if (order == 4) then
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Calculate damping flux
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          gy(i,j,k) = -1 / mesh%full_cos_lat(j)**2 * (f(i,j,k) - f(i,j-1,k)) / mesh%dlat(j)  -    &
                      mesh%full_sin_lat(j) / mesh%full_cos_lat(j) *                               &
                      (f(i,j+1,k) - f(i,j,k) - f(i,j-1,k) + f(i,j-2,k)) / (2 * mesh%dlat(j)**2) + &
                      (-f(i,j-2,k) + 3 * f(i,j-1,k) - 3 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**3
        end do
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j,k) - f(i,j-1,k))))
        end do
      end do
      call fill_halo(block, gy, full_lon=.false., full_lat=.true., full_lev=.true.)
      ! Damp physical variable at last.
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          f(i,j,k) = f(i,j,k) - dt * (gy(i,j+1,k) - gy(i,j,k)) / mesh%dlat(j) * cy_half_lat(j,k,4)
        end do
      end do
    end do
    else if (order == 2) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
!            f(i,j,k) = f(i,j,k) + dt * (-mesh%half_sin_lat(j) / mesh%half_cos_lat(j) * (f(i,j+1,k) - f(i,j-1,k)) /&
!                      (2.d0 * mesh%dlat(j)) +(f(i,j-1,k) - 2.d0 * f(i,j,k) + f(i,j+1,k)) / mesh%dlat(j)**2) * cy_half_lat(j,k,2)
            f(i,j,k) = f(i,j,k) + dt / mesh%half_cos_lat(j) / mesh%dlat(j)**2 * &
                        (mesh%full_cos_lat(j+1) * (f(i,j+1,k) - f(i,j  ,k)) - &
                         mesh%full_cos_lat(j  ) * (f(i,j  ,k) - f(i,j-1,k))) * cy_half_lat(j,k,2)
          end do
        end do
      end do
    end if

    call fill_halo(block, f, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine meridional_damp_on_vtx

end module meridional_damp_mod
