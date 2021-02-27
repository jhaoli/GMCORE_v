module vor_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use pv_mod
  use reduce_mod

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp_run

  real(r8), allocatable :: cv_full_lat(:,:)
  real(r8), allocatable :: cv_half_lat(:,:)

contains

  subroutine vor_damp_init()

    integer j, jr0, jr, k

    call vor_damp_final()

    ! Only do vorticity damping in reduced regions.
    ! First, find the interface when reduce starts.
    jr0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        if (reduce_factors(jr) > 1) jr0 = jr
      end if
    end do

    allocate(cv_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cv_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (vor_damp_order)
    case (2)
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          if (global_mesh%full_lat(j) <= 0) then
            jr = j - global_mesh%full_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%full_lat_iend_no_pole - j + 1
          end if
          cv_full_lat(j,k) = vor_damp_coef2 * exp(jr**2 * log(0.8_r8) / jr0**2) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          if (global_mesh%half_lat(j) <= 0) then
            jr = j - global_mesh%half_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%half_lat_iend_no_pole - j + 1
          end if
          cv_half_lat(j,k) = vor_damp_coef2 * exp(jr**2 * log(0.8_r8) / jr0**2) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case default
      call log_error('Unsupported vor_damp_order ' // trim(to_str(vor_damp_order)) // '!')
    end select

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(cv_full_lat)) deallocate(cv_full_lat)
    if (allocated(cv_half_lat)) deallocate(cv_half_lat)

  end subroutine vor_damp_final

  subroutine vor_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    select case (vor_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            state%u(i,j,k) = state%u(i,j,k) - dt * cv_full_lat(j,k) * ( &
              state%vor(i,j+1,k) - state%vor(i,j,k)) / mesh%le_lon(j)
#else
            state%u(i,j,k) = state%u(i,j,k) - dt * cv_full_lat(j,k) * ( &
              state%vor(i,j,k) - state%vor(i,j-1,k)) / mesh%le_lon(j)
#endif
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
#else
            state%v(i,j,k) = state%v(i,j,k) + dt * cv_half_lat(j,k) * ( &
              state%vor(i,j,k) - state%vor(i-1,j,k)) / mesh%le_lat(j)
#endif
          end do
        end do
      end do
    case (4)
    end select

  end subroutine vor_damp_run

end module vor_damp_mod
