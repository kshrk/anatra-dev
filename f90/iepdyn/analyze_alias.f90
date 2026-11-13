!-------------------------------------------------------------------------------
    subroutine add_aliased_states(option, boundary, f, fwrk)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option),   intent(inout) :: option
      type(s_boundary), intent(inout) :: boundary
      type(s_func),     intent(inout) :: f
      type(s_fwrk),     intent(inout) :: fwrk

      ! Local
      !
      integer :: nstate_old, nstate, nalias
      integer :: i, j, k, ik, jk
      integer :: ialias
      integer :: nt_range, nboundary

      ! Arrays
      !
      logical, allocatable :: is_connected_old(:, :)
      real(8), allocatable :: h_old(:, :)
      real(8), allocatable :: R_old(:, :, :)
      integer, allocatable :: alias_of(:)


      ! Setup
      !
      nstate_old = option%nstate
      nalias     = option%nalias
      nstate     = nstate_old + nalias
      nt_range   = option%nt_range 

      ! Allocate
      !
      allocate(is_connected_old(nstate_old, nstate_old))
      allocate(h_old(nstate_old, nstate_old))
      allocate(R_old(0:nt_range, nstate_old, nstate_old))
      allocate(alias_of(nstate))

      ! Backup
      !
      is_connected_old = boundary%is_connected
      h_old            = fwrk%h
      R_old            = f%R

      ! Define mapping from new states to original states
      !
      do i = 1, nstate_old
        alias_of(i) = i
      end do

      do ialias = 1, nalias
        i = nstate_old + ialias
        alias_of(i) = option%alias_state_ids(ialias)
      end do

      ! Append is_connected
      !
      deallocate(boundary%is_connected)
      allocate(boundary%is_connected(nstate, nstate))

      do i = 1, nstate
        do j = 1, nstate
          boundary%is_connected(i, j) = &
            is_connected_old(alias_of(i), alias_of(j))
        end do
      end do

      ! Update nboundary and boundary information
      !
      option%nstate = nstate
      call define_boundary(option, boundary)

      ! Append h
      !
      deallocate(fwrk%h)
      allocate(fwrk%h(nstate, nstate))

      do i = 1, nstate
        do j = 1, nstate
          fwrk%h(i, j) = h_old(alias_of(i), alias_of(j))
        end do
      end do

      ! Append kmesh and K
      !
      do i = 1, nstate
        do j = 1, nstate
          do k = 1, nstate

            ! Only append newly introduced combinations
            !
            if (i <= nstate_old .and. &
                j <= nstate_old .and. &
                k <= nstate_old) cycle

            ik = fwrk%kmesh(alias_of(i), alias_of(j), alias_of(k))

            if (ik /= 0) then
              fwrk%nk  = fwrk%nk + 1
              jk       = fwrk%nk

              fwrk%kmesh(i, j, k) = jk
              fwrk%K(:, jk)       = fwrk%K(:, ik)
            end if

          end do
        end do
      end do

      ! Rebuild K, M, and hit_count
      !
      deallocate(f%K, f%M, f%hit_count)
      call convert_Kijk_arrays(option, boundary, fwrk, f)
      call normalize_Kijk(option, boundary, f, verbose = .false.)
      call calc_Mjk_from_Kijk(option, boundary, f)

      ! Append R
      !
      deallocate(f%R, f%Rint, f%P0)
      allocate(f%R   (0:nt_range, nstate, nstate))
      allocate(f%Rint(0:nt_range, nstate, nstate))
      allocate(f%P0  (0:nt_range, nstate))
      f%R    = 0.0d0
      f%Rint = 0.0d0
      f%P0   = 0.0d0

      do i = 1, nstate
        do j = 1, nstate
          f%R(:, i, j) = R_old(:, alias_of(i), alias_of(j))
        end do
      end do
      call running_integral_Rij(option, f)
      call calc_P0_from_Rij    (option, boundary, f)
       

      ! Deallocate
      !
      deallocate(is_connected_old, h_old, R_old)
      deallocate(alias_of)
!
    end subroutine add_aliased_states
!-------------------------------------------------------------------------------
