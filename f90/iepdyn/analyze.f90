!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl

  ! Constants
  !

  ! Structures
  !
  type :: s_state
    integer :: nmol
    integer :: nstep
    integer :: unperturbed_id = -1
    integer :: use_for_Rij    = .true.
    integer, allocatable :: data(:, :) ! (nstep, nmol)
  end type s_state

  type :: s_boundary
    integer :: nboundary
    logical, allocatable :: is_connected(:, :)
    integer, allocatable :: p2b(:, :)  ! pair        => boundary id
    integer, allocatable :: b2p(:, :)  ! boundary id => pair
    integer, allocatable :: n_influx_boundary(:)
    integer, allocatable :: influx_boundary(:, :)
    logical, allocatable :: conv_direc(:)
  end type s_boundary

  type :: s_func
    real(8), allocatable :: K(:, :, :)
    real(8), allocatable :: M(:, :)
    real(8), allocatable :: P0(:, :)
    real(8), allocatable :: R(:, :, :), Rint(:, :, :)
    real(8), allocatable :: hit_count(:) 
  end type s_func

  type :: s_inpcond
    integer, allocatable :: nfile_each_state(:)
    integer, allocatable :: nfile_to_be_read(:)
  end type s_inpcond

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, einput, output, option, timegrid)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),       intent(in)    :: input
      type(s_extra_input), intent(in)    :: einput
      type(s_output),      intent(in)    :: output
      type(s_option),      intent(in)    :: option
      type(s_timegrid),    intent(in)    :: timegrid 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname 

      ! Local
      !
      character(len=MaxChar) :: line
      integer                :: ndim, nstep, nfile, nstate
      integer                :: nt_range
      integer                :: nboundary
      logical                :: is_end

      type(s_func)     :: f
      type(s_inpcond)  :: ic
      type(s_boundary) :: boundary
      type(s_cv)       :: cv
      type(s_state)    :: state

      ! Dummy
      !
      integer :: ifile, istep, iseg, is, js, is1, is2, ib, id, idir

      ! Arrays
      !
      integer, allocatable :: unperturbed_ids(:), use_for_Rij(:)
      real(8), allocatable :: Ktmp(:, :, :, :), htmp(:, :)


      ! Setup
      !
      ndim     = option%ndim * option%nmol 
      nfile    = input%ncv
      nstate   = option%nstate
      nt_range = option%nt_range

      ! Allocate 
      !
      allocate(f%R   (0:nt_range, nstate, nstate))
      allocate(f%Rint(0:nt_range, nstate, nstate))
      allocate(f%P0  (0:nt_range, nstate))

      allocate(Ktmp(0:nt_range, nstate, nstate, nstate))
      allocate(htmp(nstate, nstate))

      f%R    = 0.0d0
      f%Rint = 0.0d0
      f%P0   = 0.0d0
      Ktmp   = 0.0d0
      htmp   = 0.0d0
      
      ! Read Unperturbed_ID file
      !
      if (option%use_perturbed_traj) then
        write(iw,*)
        write(iw,'("Analyze> Read f_unperturbed_id file")')
        write(iw,'("Note: unperturbed state info. is used only if use_perturbed_traj = .true.")')

        ! Get unperturbed state id for each file
        !
        allocate(unperturbed_ids(nfile), use_for_Rij(nfile))
        call read_f_unperturbed_id(option, nfile, unperturbed_ids, use_for_Rij)
        do ifile = 1, nfile
          write(iw,'(3i10)') ifile, unperturbed_ids(ifile), use_for_Rij(ifile)
        end do
       
        ! Get # of files for each unperturbed state
        !
        allocate(ic%nfile_each_state(nstate))
        allocate(ic%nfile_to_be_read(nstate))

        ic%nfile_each_state = 0
        do ifile = 1, nfile
          is                      = unperturbed_ids(ifile)
          ic%nfile_each_state(is) = ic%nfile_each_state(is) + 1
        end do

        ! Initial setting (will not be changed if check_senserr = .false.)
        !
        do is = 1, nstate
          ic%nfile_to_be_read(is) = ic%nfile_each_state(is)  
        end do

      end if

      ! TODO: Following should be capsuled as a subroutine in future update 
      !
      if (option%input_type == InputTypeTIMESERIES) then

        do ifile = 1, nfile
       
          write(iw,'("Analyze> Read CV file: ", 2x,a)') trim(input%fcv(ifile))
       
          iseg   = 0
          is_end = .false.
          call open_file(input%fcv(ifile), io, stat = 'old')
          do while (.not. is_end)
       
            ! Initialized
            ! 
            if (allocated(cv%x)) then
              call deallocate_cv(cv)
            end if
       
            if (allocated(state%data)) then
              state%nmol  = 0
              state%nstep = 0
              deallocate(state%data)
            end if
       
            ! Read
            !
            call read_cv_split(io, ndim, '#SPLIT', cv, is_end)
       
            if (cv%nstep /= 0) then
              iseg = iseg + 1
              call get_state(option, cv, state)
       
              state%unperturbed_id = -1
              state%use_for_Rij    = -1
              if (option%use_perturbed_traj) then
                state%unperturbed_id = unperturbed_ids(ifile)
                state%use_for_Rij    = use_for_Rij(ifile)
              end if
       
              ! Update Connectivity
              !
              call get_state_connectivity(output, option, state, boundary)
       
              ! Update R- and K-functions
              !
              call update_Rij_wo_normalize(option, state, f)
              call update_Kijk_wo_normalize(option, state, Ktmp, htmp)
            end if
       
          end do
       
          close(io)
       
        end do

        if (option%output_histogram) then
          call output_Rij_hist(option, output, f)
          call output_Kijk_hist(option, output, Ktmp, htmp)
        end if

      else if (option%input_type == InputTypeHISTOGRAM) then
        do ifile = 1, nfile
          write(iw,'("Analyze> Read CV file: ", 2x,a)') trim(input%fcv(ifile))
          iseg   = 0
          is_end = .false.
          call open_file(input%fcv(ifile), io, stat = 'old')
          read(io,'(a)') line
          if (trim(line) == 'KIJK') then
            call update_Kijk_from_hist(io, option, ktmp, htmp) 
          else if (trim(line) == 'RIJ') then
            call update_Rij_from_hist(io, option, f) 
          end if
          close(io)
        end do
        call get_state_connectivity_from_h(option, htmp, boundary)
      end if

      ! Show connectivity
      !
      write(iw,*)
      write(iw,'("Analyze> Get Connectivity")')
      call show_state_connectivity(option, boundary)

      ! Define boundary
      !
      write(iw,*)
      write(iw,'("Analyze> Define Boundary")')
      call define_boundary(option, boundary)
      nboundary = boundary%nboundary
      write(iw,'(">> Done")')

      ! Convert arrays
      !
      allocate(f%K(0:nt_range, nstate, -nboundary:nboundary))
      allocate(f%M(0:nt_range, -nboundary:nboundary))
      allocate(f%hit_count(-nboundary:nboundary))

      f%K         = 0.0d0
      f%M         = 0.0d0
      f%hit_count = 0.0d0

      call convert_Kijk_arrays(option, boundary, Ktmp, htmp, f)

      ! Normalize R- and K-functions
      ! 
      call normalize_Rij(option, f)
      call normalize_Kijk(option, boundary, f)

      ! Check Kijk
      !
      call check_Kijk(option, boundary, f)

      ! Compute R-integration and P0
      !
      write(iw,*)
      write(iw,'("Analyze> Calculate P0 from Rij")')

      call running_integral_Rij(option, f)
      call calc_P0_from_Rij    (option, boundary, f)
      write(iw,'(">> Done")')

      ! Compute Mij
      !
      write(iw,*)
      write(iw,'("Analyze> Calculate Mjk from Kijk")')
      call calc_Mjk_from_Kijk(option, boundary, f)
      write(iw,'(">> Done")')

      ! Output
      !
      write(iw,*)
      write(iw,'("Analyze> Print out TCFs")')
      call write_Rij   (output, option, boundary, f)
      call write_P0    (output, option, f)
      call write_Kijk  (output, option, boundary, f)
      call write_Mjk   (output, option, boundary, f)
      write(iw,'(">> Done")')

      ! Setup Boundary conditions 
      !
      write(iw,*)
      write(iw,'("Analyze> Setup boundary conditions")')
      call set_reflection(output, option, boundary, f) 
      call set_product   (output, option, boundary, f)
      write(iw,'(">> Done")')

      if (option%calc_Pint .or. option%calc_Steady) then
        call reacdyn_pint(output, option, boundary, f)
      end if

      if (option%extrapolate) then

        write(iw,*)
        write(iw,'("Analyze> Start propagation")')

        ! Extend timescale 
        !
        call reacdyn_tcf(output, option, boundary, f)

        write(iw,'(">> Done")')

      end if 

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(option, cv, state) 
!-----------------------------------------------------------------------
      implicit none
      
      type(s_option) :: option
      type(s_cv)     :: cv
      type(s_state)  :: state

      ! Local
      !
      integer :: nmol, ndim, nstep 

      ! Dummy
      !
      integer :: istep, imol, idim, jdim, ia, istate

      ! Array
      !
      real(8), allocatable :: val(:)


      ! Setup
      !
      nmol  = option%nmol
      ndim  = option%ndim
      nstep = cv%nstep

      state%nmol  = nmol
      state%nstep = nstep

      allocate(state%data(nstep, nmol))
      allocate(val(ndim))

      state%data = - 1 
      do istep = 1, nstep
        jdim = 0
        do imol = 1, nmol
          do idim = 1, ndim
            jdim      = jdim + 1
            val(idim) = cv%data(jdim, istep)
          end do

          do istate = 1, option%nstate

            ia = 0
            do idim = 1, ndim
              if (val(idim) >= option%state_def(1, idim, istate) &
            .and. val(idim) <  option%state_def(2, idim, istate)) then
                ia = ia + 1
              end if 
            end do

            if (ia == ndim) then
              state%data(istep, imol) = istate  
            end if

          end do

          if (state%data(istep, imol) == - 1) then
            write(iw,'("Get_State> Error.")')
            write(iw,'("State determination failed. &
                       &Please check state definition.")')
            stop
          end if

        end do
      end do 

    end subroutine get_state
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_f_unperturbed_id(option, nfile, unperturbed_ids, use_for_Rij)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: nfile 
      integer,        intent(inout) :: unperturbed_ids(nfile)
      integer,        intent(inout) :: use_for_Rij(nfile) 

      ! I/O
      !
      integer :: io

      ! Dummy
      !
      integer :: ifile


      ! Initialize 
      !
      unperturbed_ids = 0
      use_for_Rij     = 0

      ! Read
      !
      call open_file(option%f_unperturbed_id, io, stat = 'old')
      do ifile = 1, nfile
        read(io,*) unperturbed_ids(ifile), use_for_Rij(ifile) 
      end do
      close(io)
      
!
    end subroutine read_f_unperturbed_id 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state_connectivity(output, option, state, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_state),    intent(in)    :: state 
      type(s_boundary), intent(inout) :: boundary 

      ! Local
      !
      integer :: nstate, nmol, nstep

      ! Dummy
      !
      integer :: istep, imol
      integer :: snow, sprev

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname


      ! Setup
      !
      nmol   = option%nmol
      nstate = option%nstate
      nstep  = state%nstep

      if (.not. allocated(boundary%is_connected)) then
        allocate(boundary%is_connected(nstate, nstate))
        boundary%is_connected = .false.
      end if

      do imol = 1, nmol
        do istep = 1, nstep
          snow = state%data(istep, imol)
          if (istep == 1) then
            sprev = snow
            cycle
          end if

          if (snow /= sprev) then
            boundary%is_connected(snow,  sprev) = .true. 
            boundary%is_connected(sprev, snow)  = .true. 
          end if

          sprev = snow

        end do
      end do

    end subroutine get_state_connectivity
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state_connectivity_from_h(option, h, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option 
      real(8),          intent(in)    :: h(option%nstate, option%nstate)
      type(s_boundary), intent(inout) :: boundary 

      ! Local
      !
      integer :: nstate

      ! Dummy
      !
      integer :: js, ks 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname


      ! Setup
      !
      nstate = option%nstate

      if (.not. allocated(boundary%is_connected)) then
        allocate(boundary%is_connected(nstate, nstate))
        boundary%is_connected = .false.
      end if

      do ks = 1, nstate
        do js = 1, nstate
          if (h(js, ks) > 0.999d0) then
            boundary%is_connected(js,  ks) = .true. 
            boundary%is_connected(ks,  js) = .true. 
          end if
        end do 
      end do

    end subroutine get_state_connectivity_from_h
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_state_connectivity(option, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary 

      ! Local
      !
      integer :: nstate, nmol

      ! Dummy
      !
      integer :: istep, imol, is, js
      integer :: snow, sprev


      ! Setup
      !
      nmol   = option%nmol
      nstate = option%nstate

      do js = 1, nstate - 1
        do is = js + 1, nstate
          if (boundary%is_connected(is, js)) then
            write(iw,'(i0,2x,i0)') is, js
          end if
        end do 
      end do

    end subroutine show_state_connectivity
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine define_boundary(option, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary 

      ! Local
      !
      integer :: nstate, nboundary

      ! Dummy
      !
      integer :: is, js, ib, ic 


      ! Setup
      !
      nstate = option%nstate

      allocate(boundary%p2b(nstate, nstate))

      boundary%p2b = 0
      nboundary    = 0
      do js = 1, nstate - 1
        do is = js + 1, nstate
          if (boundary%is_connected(is, js)) then
            nboundary            =  nboundary + 1
            boundary%p2b(is, js) =  nboundary
            boundary%p2b(js, is) = -nboundary
          end if
        end do
      end do

      boundary%nboundary = nboundary

      allocate(boundary%b2p(2, -nboundary:nboundary))

      nboundary    = 0
      boundary%b2p = 0
      do js = 1, nstate - 1
        do is = js + 1, nstate
          if (boundary%is_connected(is, js)) then
            nboundary                   = nboundary + 1
            boundary%b2p(1,  nboundary) = js
            boundary%b2p(2,  nboundary) = is
            boundary%b2p(1, -nboundary) = is
            boundary%b2p(2, -nboundary) = js
          end if
        end do
      end do

      ! Setup influx_boundary
      !
      nboundary = boundary%nboundary

      allocate(boundary%n_influx_boundary(nstate))
      allocate(boundary%influx_boundary(2*nboundary, nstate))

      do is = 1, nstate
        ic = 0
        do ib = -nboundary, nboundary
          if (ib == 0) cycle
          js = boundary%b2p(2, ib)
          if (js == is) then
            ic = ic + 1
            boundary%influx_boundary(ic, is) = ib 
          end if
        end do
        boundary%n_influx_boundary(is) = ic
      end do

    end subroutine define_boundary 
!-----------------------------------------------------------------------

    include 'analyze_Rij.f90'
    include 'analyze_P0.f90'
    include 'analyze_Kijk.f90'
    include 'analyze_Mij.f90'
    include 'analyze_bc.f90'
    include 'analyze_reacdyn.f90'
    include 'analyze_pint.f90'

end module mod_analyze
!=======================================================================
