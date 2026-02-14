!-----------------------------------------------------------------------
    subroutine update_Kijk_wo_normalize(option, state, Kijk, hit_count)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_state),    intent(in)    :: state
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             option%nstate,             &
                                             option%nstate)
      real(8),          intent(inout) :: hit_count(option%nstate,       &
                                                   option%nstate)

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse
      integer :: init_id, unp_id
      real(8) :: dt
      logical :: is_final, is_prod, is_dissoc, use_single_event

      ! Dummy 
      !
      integer :: istep, jstep, imol
      integer :: is, js, ks, ls, ib, jb
      integer :: it, it_reac, it_diff, nt
      integer :: ista, iend, ireac
      integer :: iprod, idissoc


      ! Setup
      !
      nmol             = option%nmol
      nstep            = state%nstep
      unp_id           = state%unperturbed_id
      nt_range         = option%nt_range
      nt_sparse        = option%nt_sparse
      dt               = option%dt_out
      use_single_event = option%use_single_event

      do imol = 1, nmol

        ! Search reaction time
        !
        ks = state%data(1, imol)
        js = ks

        if (allocated(state%init_id)) then
          init_id  = state%init_id(imol) 
        else
          init_id  = ks
        end if

        ireac    = 0
        is_final = .false. 
        do istep = nt_sparse + 1, nstep, nt_sparse
          is = state%data(istep, imol)
          if (is /= js) then
            ireac  = ireac + 1
            
            if (ireac == 1) then
              ista              = istep - nt_sparse
              js                = is
              if (unp_id == js .or. unp_id == -1) then 
                hit_count(js, ks) = hit_count(js, ks) + 1.0d0
              end if
            else
              iend    = istep - nt_sparse 
              it_diff = iend - ista
              it_diff = int(it_diff / dble(nt_sparse))

              if (unp_id == js .or. unp_id == -1) then 
                Kijk(it_diff, is, js, ks) &
                        = Kijk(it_diff, is, js, ks) + 1.0d0
              end if

              if (unp_id == is .or. unp_id == -1) then
                hit_count(is, js) = hit_count(is, js) + 1.0d0
              end if

              ls      = ks
              ks      = js
              js      = is
              ista    = istep - nt_sparse 

              !if (istep == nstep) then
              if (istep == nstep .or. istep + nt_sparse > nstep) then
                is_final = .true.
              end if

              if (use_single_event) then
                if (ls /= init_id) then
                  Kijk(it_diff, js, ks, ls) &
                    = Kijk(it_diff, js, ks, ls)    - 1.0d0
                  hit_count(ks, ls) = hit_count(ks, ls) - 1.0d0
                end if
              end if

            end if

          end if
        end do

        ! New implementation
        !
        if (ireac == 0) cycle

        if (option%use_product_state) then
          is_prod = .false. 
          do iprod = 1, option%nproduct
            if (js == option%product_state_ids(iprod)) then
              is_prod = .true.
              exit 
            end if
          end do
          if (is_prod) cycle
        end if
!
!        if (option%use_dissociate_state) then
!          is_dissoc = .false.
!          do idissoc = 1, option%ndissoc
!            if (js == option%dissociate_state_ids(idissoc)) then
!              is_dissoc = .true.
!              exit 
!            end if
!          end do
!          if (is_dissoc) cycle 
!        end if
!
!        if (.not. option%is_dissoc(js)) then
!          if (unp_id == js .or. unp_id == -1) then
!            hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!          end if
!        end if

        if (option%is_dissoc(js)) then
          if (is_final .and. (unp_id == js .or. unp_id == -1)) then
            hit_count(js, ks) = hit_count(js, ks) - 1.0d0
          end if
        else
          if (unp_id == js .or. unp_id == -1) then
            hit_count(js, ks) = hit_count(js, ks) - 1.0d0
          end if
        end if

        ! End of New implementation
        !

!        if (ireac == 0 .or. ireac == 1) then
!
!          if (option%use_product_state) then
!            do iprod = 1, option%nproduct
!              if (js == option%product_state_ids(iprod)) then
!                return
!              end if
!            end do 
!          end if
!
!          if (option%use_dissociate_state) then
!            do idissoc = 1, option%ndissoc
!              if (js == option%dissociate_state_ids(idissoc)) then
!                return
!              end if
!            end do 
!          end if
!
!          if (ireac == 1) then
!            if (.not. option%is_dissoc(js)) then
!              if (unp_id == js .or. unp_id == -1) then
!                hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!              end if
!            end if
!            return
!          end if
!
!          !write(iw,'("Update_Kijk_wo_normalize> Error.")')
!          !write(iw,'("No reaction is observed.")')
!          !stop
!
!          !if (option%use_product_state) then
!          !  do iprod = 1, option%nproduct
!          !    if (js == option%product_state_ids(iprod)) then
!          !      return
!          !    end if
!          !  end do 
!          !else if (option%use_dissociate_state) then
!          !  do idissoc = 1, option%ndissoc
!          !    if (js == option%dissociate_state_ids(idissoc)) then
!          !      return
!          !    end if
!          !  end do 
!          !else
!          !  write(iw,'("Update_Kijk_wo_normalize> Error.")')
!          !  write(iw,'("No reaction is observed.")')
!          !  stop
!          !end if
!
!          !if (ireac == 1) then
!          !  if (.not. option%is_dissoc(js)) then 
!          !    hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!          !  end if
!          !end if
!
!          !if (ireac == 1) then
!          !  if (.not. option%is_dissoc(js)) then
!          !    if (unp_id == js .or. unp_id == -1) then
!          !      hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!          !    end if
!          !  end if
!          !end if
!
!        else
!
!          if (.not. option%is_dissoc(js)) then
!            if (unp_id == js .or. unp_id == -1) then
!              hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!            end if
!          end if
!
!        end if

      end do

    end subroutine update_Kijk_wo_normalize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine convert_Kijk_arrays(option, boundary, Ktmp, htmp, Kijk, hit_count)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      real(8),          intent(in)    :: Ktmp(0:option%nt_range, &
                                              option%nstate,     &
                                              option%nstate,     &
                                              option%nstate)
      real(8),          intent(in)    :: htmp(option%nstate, option%nstate)
      real(8),          intent(inout) :: Kijk(0:option%nt_range, &
                                              option%nstate,     &
                                              -boundary%nboundary:boundary%nboundary)
      real(8),          intent(inout) :: hit_count(-boundary%nboundary:boundary%nboundary)

      ! Local
      integer :: nboundary

      ! Dummy
      integer :: ib, jb, is1, is2, js1, js2
      integer :: inflx


      nboundary = boundary%nboundary
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        is1           = boundary%b2p(1, ib)
        is2           = boundary%b2p(2, ib)
        hit_count(ib) = htmp(is2, is1)

        do inflx = 1, boundary%n_influx_boundary(is1)
          jb  = boundary%influx_boundary(inflx, is1)
          js1 = boundary%b2p(1, jb)
          js2 = boundary%b2p(2, jb)
          Kijk(:, is2, jb) = Ktmp(:, is2, is1, js1) 
        end do 
      end do

!
    end subroutine convert_Kijk_arrays
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_Kijk(option, boundary, Kijk, hit_count)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)
      real(8),          intent(inout) :: hit_count(-boundary%nboundary: &
                                                    boundary%nboundary)
                                           


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, is1, is2, js, ib
      real(8) :: state_sum


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      dt        = option%dt_out

      ! Normalize
      !
      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        state_sum      = hit_count(ib)

        if (state_sum < 1.0d-10) then
          Kijk(:, :, ib) = 0.0d0
        else 
          Kijk(:, :, ib) = Kijk(:, :, ib) / (state_sum * dt) 
        end if

      end do

      ! Report Hit counts
      !
      write(iw,'("Normalize_Kijk> Summary of Hit counts")')
      write(iw,'("I <--- J")')
      do ib = -nboundary, nboundary
        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        write(iw,'(i5,i5," : ", f20.10)') is2, is1, hit_count(ib) 
      end do
      write(iw,*)
      

    end subroutine normalize_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine renormalize_Kijk(option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, is1, is2, js, ib
      real(8) :: state_sum


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      dt        = option%dt_out

      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        if (.not. option%is_dissoc(is2)) then 
          state_sum = 0.0d0
          do js = 1, nstate
            if (.not. boundary%is_connected(js, is2)) cycle 
            do istep = 0, nt_range
              state_sum = state_sum + Kijk(istep, js, ib)
            end do
          end do

          if (state_sum > 1.0d-10) then
            Kijk(:, :, ib) = Kijk(:, :, ib) / (state_sum * dt)
          else
            Kijk(:, :, ib) = 0.0d0
          end if 

        end if

      end do

    end subroutine renormalize_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Kijk(output, option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      real(8),          intent(in) :: Kijk(0:option%nt_range,  &
                                           option%nstate,      &
                                          -boundary%nboundary:boundary%nboundary)

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range, nboundary
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer :: is, js, is1, is2, ib, id, istep


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nboundary = boundary%nboundary


      write(fname, '(a,".Kijk")') trim(output%fhead)
      call open_file(fname, io)
      
      write(io,'("# J <--- B")')
      id = 1 
      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        do js = 1, nstate
          if (boundary%is_connected(js, is2)) then
            write(io,'(" # Col. ", i0, " Kijk : ", i0, " <= ", i0, " < ", i0)')  &
              id + 1, js, is2, is1

            id = id + 1
          end if

        end do

      end do

      do istep = 0, nt_range
        write(io,'(f20.10)', advance = 'no') option%dt_out * istep 
        do ib = -nboundary, nboundary

          if (ib == 0) then
            cycle
          end if

          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)

          do js = 1, nstate
            if (boundary%is_connected(js, is2)) then
              write(io,'(e15.7,2x)', advance = 'no') Kijk(istep, js, ib)
            end if 
          end do
        end do
        write(io,*)
      end do
      close(io)

    end subroutine write_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_Kijk(option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)
                                           


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, is1, is2, js, ib
      real(8) :: state_sum


      if (.not. option%check_Kijk) return

      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      dt        = option%dt_out

      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        state_sum = 0.0d0
        do js = 1, nstate
          if (boundary%is_connected(js, is2)) then
            state_sum = state_sum + sum(Kijk(:, js, ib)) * dt
            !if (state_sum < 1.0d-6) then
            !  write(iw,'("Check_Kijk> Warning.")')
            !  write(iw,'("K(",i0,2x,i0,2x,i0,") is zero")') js, is2, is1
            !end if 
          end if
        end do

        write(iw,'("Kint value at Boundary ", i5, i5, " : ", f15.7)') is2, is1, state_sum

      end do

    end subroutine check_Kijk
!-----------------------------------------------------------------------
