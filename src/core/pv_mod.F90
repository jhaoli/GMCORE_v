module pv_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public diag_pv
  public interp_pv_init
  public interp_pv_final
  public interp_pv_midpoint
  public interp_pv_upwind
  public interp_pv_apvm

  real(r8), allocatable :: upwind_wgt_lon(:), &
                           upwind_wgt_lat(:)

contains
  
  subroutine interp_pv_init()

    integer j
    
    call interp_pv_final()

    allocate(upwind_wgt_lon(global_mesh%num_full_lat))
    allocate(upwind_wgt_lat(global_mesh%num_half_lat))

    do j = 1, global_mesh%num_full_lat
      if (abs(global_mesh%full_lat_deg(j)) > 80) then
        upwind_wgt_lon(j) = 1 / upwind_wgt_pv
      else
        upwind_wgt_lon(j) = 1 + (1 / upwind_wgt_pv - 1) * &
          exp(-0.1 * (abs(global_mesh%full_lat_deg(j)) - 80)**2)
      end if
    end do

    do j = 1, global_mesh%num_half_lat
      if (abs(global_mesh%half_lat_deg(j)) > 80) then
        upwind_wgt_lat(j) = 1 / upwind_wgt_pv
      else
        upwind_wgt_lat(j) = 1 + (1 / upwind_wgt_pv - 1) * &
          exp(-0.1 * (abs(global_mesh%half_lat_deg(j)) - 80)**2)
      end if
    end do

  end subroutine interp_pv_init
  
  subroutine interp_pv_final()

    if (allocated(upwind_wgt_lon)) deallocate(upwind_wgt_lon)
    if (allocated(upwind_wgt_lat)) deallocate(upwind_wgt_lat)

  end subroutine interp_pv_final

  subroutine calc_vor(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole(state%mesh%num_full_lev)
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%vor(i,j,k) = (                                                                &
            state%u(i  ,j-1,k) * mesh%de_lon(j-1) - state%u(i  ,j  ,k) * mesh%de_lon(j  ) + &
            state%v(i+1,j  ,k) * mesh%de_lat(j  ) - state%v(i  ,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
#else
          state%vor(i,j,k) = (                                                                 &
            state%u(i  ,j  ,k) * mesh%de_lon(j  ) - state%u(i  ,j+1,k) * mesh%de_lon(j+1) + &
            state%v(i+1,j  ,k) * mesh%de_lat(j  ) - state%v(i  ,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
#endif
        end do
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole(k) = pole(k) - state%u(i,j,k) * mesh%de_lon(j)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_half_lon / mesh%area_vtx(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%vor(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole(k) = pole(k) + state%u(i,j-1,k) * mesh%de_lon(j-1)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_half_lon / mesh%area_vtx(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%vor(i,j,k) = pole(k)
        end do
      end do
    end if
#else
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            pole(k) = pole(k) - state%u(i,j+1,k) * mesh%de_lon(j+1)
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%vor(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            pole(k) = pole(k) + state%u(i,j,k) * mesh%de_lon(j)
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%vor(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
#endif
    call fill_halo(block, state%vor, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine calc_vor

  subroutine diag_pv(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole(state%mesh%num_full_lev)
    integer i, j, k

    call calc_vor(block, state)

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%pv(i,j,k) = (state%vor(i,j,k) + mesh%half_f(j)) / state%m_vtx(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, state%pv, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine diag_pv

  subroutine calc_dpv_edge(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    ! Tangent pv difference
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%dpv_lat_t(i,j,k) = state%pv(i,j,k) - state%pv(i-1,j,k)
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%dpv_lon_t(i,j,k) = state%pv(i,j+1,k) - state%pv(i,j  ,k)
#else
          state%dpv_lon_t(i,j,k) = state%pv(i,j  ,k) - state%pv(i,j-1,k)
#endif
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

    ! Normal pv difference
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%dpv_lat_n(i,j,k) = 0.25_r8 * (state%dpv_lon_t(i-1,j-1,k) + state%dpv_lon_t(i,j-1,k) + &
                                              state%dpv_lon_t(i-1,j  ,k) + state%dpv_lon_t(i,j  ,k))
#else
          state%dpv_lat_n(i,j,k) = 0.25_r8 * (state%dpv_lon_t(i-1,j  ,k) + state%dpv_lon_t(i,j  ,k) + &
                                              state%dpv_lon_t(i-1,j+1,k) + state%dpv_lon_t(i,j+1,k))
#endif
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%dpv_lon_n(i,j,k) = 0.25_r8 * (state%dpv_lat_t(i,j  ,k) + state%dpv_lat_t(i+1,j  ,k) + &
                                              state%dpv_lat_t(i,j+1,k) + state%dpv_lat_t(i+1,j+1,k))
#else
          state%dpv_lon_n(i,j,k) = 0.25_r8 * (state%dpv_lat_t(i,j-1,k) + state%dpv_lat_t(i+1,j-1,k) + &
                                              state%dpv_lat_t(i,j  ,k) + state%dpv_lat_t(i+1,j  ,k))
#endif
        end do
      end do
    end do

  end subroutine calc_dpv_edge

  subroutine interp_pv_midpoint(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%pv_lat(i,j,k) = 0.5_r8 * (state%pv(i-1,j,k) + state%pv(i,j,k))
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%pv_lon(i,j,k) = 0.5_r8 * (state%pv(i,j,k) + state%pv(i,j+1,k))
#else
          state%pv_lon(i,j,k) = 0.5_r8 * (state%pv(i,j,k) + state%pv(i,j-1,k))
#endif
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

  end subroutine interp_pv_midpoint

  subroutine interp_pv_upwind(block, state, upwind_wgt_, enhance_pole)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in), optional :: upwind_wgt_
    logical, intent(in), optional :: enhance_pole

    type(mesh_type), pointer :: mesh
    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8
    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8
    real(r8) beta_lon(state%mesh%full_lat_ibeg:state%mesh%full_lat_iend), &
             beta_lat(state%mesh%half_lat_ibeg:state%mesh%half_lat_iend)
    integer i, j, k
    mesh => state%mesh

    beta_lon = merge(upwind_wgt_, upwind_wgt_pv, present(upwind_wgt_))
    beta_lat = merge(upwind_wgt_, upwind_wgt_pv, present(upwind_wgt_))
    if (present(enhance_pole)) then
      if (enhance_pole) then
        beta_lon = beta_lon * upwind_wgt_lon(mesh%full_lat_ibeg:mesh%full_lat_iend)
        beta_lat = beta_lat * upwind_wgt_lat(mesh%half_lat_ibeg:mesh%half_lat_iend)
      end if
    end if

    select case (upwind_order_pv)
    case (1)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%pv_lat(i,j,k) = c11 * (state%pv(i,j,k) + state%pv(i-1,j,k)) + &
                                  c12 * (state%pv(i,j,k) - state%pv(i-1,j,k)) * &
                          beta_lat(j) * sign(1.0_r8, state%mf_lat_t(i,j,k))
          end do
        end do
      end do
    case (3)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%pv_lat(i,j,k) = c31 * (state%pv(i  ,j,k) + state%pv(i-1,j,k))  + &
                                  c32 * (state%pv(i+1,j,k) + state%pv(i-2,j,k))  + &
                                  c33 * (state%pv(i+1,j,k) - state%pv(i-2,j,k)   - &
                               3.0_r8 * (state%pv(i  ,j,k) - state%pv(i-1,j,k))) * &
                          beta_lat(j) * sign(1.0_r8, state%mf_lat_t(i,j,k))
          end do
        end do
      end do
    end select

#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    select case (upwind_order_pv)
    case (1)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            state%pv_lon(i,j,k) = c11 * (state%pv(i,j+1,k) + state%pv(i,j,k)) + &
                                  c12 * (state%pv(i,j+1,k) - state%pv(i,j,k)) * &
                          beta_lon(j) * sign(1.0_r8, state%mf_lon_t(i,j,k))
#else
            state%pv_lon(i,j,k) = c11 * (state%pv(i,j,k) + state%pv(i,j-1,k)) + &
                                  c12 * (state%pv(i,j,k) - state%pv(i,j-1,k)) * &
                          beta_lon(j) * sign(1.0_r8, state%mf_lon_t(i,j,k))
#endif
          end do
        end do
      end do
    case (3)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (mesh%is_full_lat_next_to_pole(j)) then
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
              state%pv_lon(i,j,k) = c11 * (state%pv(i,j+1,k) + state%pv(i,j,k)) + &
                                    c12 * (state%pv(i,j+1,k) - state%pv(i,j,k)) * &
                            beta_lon(j) * sign(1.0_r8, state%mf_lon_t(i,j,k))
#else
              state%pv_lon(i,j,k) = c11 * (state%pv(i,j,k) + state%pv(i,j-1,k)) + &
                                    c12 * (state%pv(i,j,k) - state%pv(i,j-1,k)) * &
                            beta_lon(j) * sign(1.0_r8, state%mf_lon_t(i,j,k))
#endif
            end do
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
              state%pv_lon(i,j,k) = c31 * (state%pv(i,j+1,k) + state%pv(i,j  ,k))  + &
                                    c32 * (state%pv(i,j+2,k) + state%pv(i,j-1,k))  + &
                                    c33 * (state%pv(i,j+2,k) - state%pv(i,j-1,k)   - &
                                 3.0_r8 * (state%pv(i,j+1,k) - state%pv(i,j  ,k))) * &
                            beta_lon(j) * sign(1.0_r8, state%mf_lon_t(i,j,k))

#else
              state%pv_lon(i,j,k) = c31 * (state%pv(i,j  ,k) + state%pv(i,j-1,k))  + &
                                    c32 * (state%pv(i,j+1,k) + state%pv(i,j-2,k))  + &
                                    c33 * (state%pv(i,j+1,k) - state%pv(i,j-2,k)   - &
                                 3.0_r8 * (state%pv(i,j  ,k) - state%pv(i,j-1,k))) * &
                            beta_lon(j) * sign(1.0_r8, state%mf_lon_t(i,j,k))
#endif
            end do
          endif
        end do
      end do
    end select
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

  end subroutine interp_pv_upwind

  subroutine interp_pv_apvm(block, state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    real(r8) u, v
    integer i, j, k

    call calc_dpv_edge(block, state)

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          u = state%mf_lat_t(i,j,k) / state%m_lat(i,j,k)
          v = state%v(i,j,k)
          state%pv_lat(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j,k) + state%pv(i-1,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lat_t(i,j,k) / mesh%le_lat(j) + &
            v * state%dpv_lat_n(i,j,k) / mesh%de_lat(j)   &
          ) * dt
        end do
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) state%pv_lat(:,mesh%half_lat_ibeg,:) = state%pv(:,mesh%half_lat_ibeg,:)
    if (mesh%has_north_pole()) state%pv_lat(:,mesh%half_lat_iend,:) = state%pv(:,mesh%half_lat_iend,:)
#endif
#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u = state%u(i,j,k)
          v = state%mf_lon_t(i,j,k) / state%m_lon(i,j,k)
#ifdef V_POLE
          state%pv_lon(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j+1,k) + state%pv(i,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lon_n(i,j,k) / mesh%de_lon(j) + &
            v * state%dpv_lon_t(i,j,k) / mesh%le_lon(j)   &
          ) * dt
#else
          state%pv_lon(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j-1,k) + state%pv(i,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lon_n(i,j,k) / mesh%de_lon(j) + &
            v * state%dpv_lon_t(i,j,k) / mesh%le_lon(j)   &
          ) * dt
#endif
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

  end subroutine interp_pv_apvm

end module pv_mod
