!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl

  ! constants
  !

  ! structures
  !
  type :: s_state
    integer :: nmol
    integer :: nstep
    integer, allocatable :: data(:, :) ! (nstep, nmol)
    integer, allocatable :: init_id(:) 
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
      integer :: ndim, nstep, nfile, nstate
      integer :: nt_range
      integer :: nboundary

      type(s_boundary) :: boundary
      type(s_cv)       :: cv
      type(s_state)    :: state

      ! Dummy
      !
      integer :: ifile, istep, is, js, is1, is2, ib, id, idir

      ! Arrays
      !
      integer,       allocatable :: init_id(:)
      real(8),       allocatable :: Rij(:, :, :), Rij_int(:, :, :)
      real(8),       allocatable :: P0(:, :)

      real(8),       allocatable :: Kijk(:, :, :), hit_count(:)
      real(8),       allocatable :: Mij(:, :)

      real(8),       allocatable :: Ktmp(:, :, :, :), htmp(:, :)


      ! Setup
      !
      ndim     = option%ndim * option%nmol 
      nfile    = input%ncv
      nstate   = option%nstate
      nt_range = option%nt_range

      ! Construct K-, R-, M-, and P0-functions
      !
      allocate(Rij(0:nt_range, nstate, nstate))
      allocate(Ktmp(0:nt_range, nstate, nstate, nstate))
      allocate(htmp(nstate, nstate))

      Rij  = 0.0d0
      Ktmp = 0.0d0
      htmp = 0.0d0
      
      ! Read init_id file
      !
      if (option%read_init_id) then 
        write(iw,*)
        write(iw,'("Analyze> Read f_init_id file")')
        write(iw,'("Note: init_id info. is used only if read_init_id = .true.")')
        allocate(init_id(option%nmol))
        call read_f_init_id(option, nfile, init_id)
      end if

      do ifile = 1, nfile

        ! Initialize 
        !
        if (allocated(cv%x)) then
          call deallocate_cv(cv)
          state%nmol  = 0
          state%nstep = 0
          deallocate(state%data)
        end if

        ! Read CV
        !
        write(iw,'("Analyze> Read CV file: ", 2x,a)') trim(input%fcv(ifile))
        call read_cv  (input%fcv(ifile), ndim, cv)
        call get_state(option, cv, state)
        if (option%read_init_id) then
          if (.not. allocated(state%init_id)) then
            allocate(state%init_id(option%nmol))
          end if
          state%init_id = init_id
        end if

        ! Update state connectivity
        !
        call get_state_connectivity(output, option, state, boundary)

        ! Update R- and K-functions
        !
        call update_Rij_wo_normalize(option, state, Rij)
        call update_Kijk_wo_normalize(option, state, Ktmp, htmp)
      end do

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
      allocate(Kijk(0:nt_range, nstate, -nboundary:nboundary))
      allocate(hit_count(-nboundary:nboundary))
      call convert_Kijk_arrays(option, boundary, Ktmp, htmp, Kijk, hit_count)

      ! Normalize R- and K-functions
      ! 
      call normalize_Rij(option, Rij)
      call normalize_Kijk(option, boundary, Kijk, hit_count)

      ! Check Kijk
      !
      call check_Kijk(option, boundary, Kijk)

      ! Compute R-integration and P0
      !
      allocate(Rij_int(0:nt_range, nstate, nstate))
      allocate(P0     (0:nt_range, nstate))
      call running_integral_Rij(option, Rij, Rij_int)
      call calc_P0_from_Rij    (option, boundary, Rij, P0)

      ! Compute Mij
      !
      allocate(Mij (0:nt_range, -nboundary:nboundary))
      call calc_Mij_from_Kijk(option, boundary, Kijk, Mij)

      ! Output
      !
      call write_Rij   (output, option, boundary, Rij, Rij_int)
      call write_P0    (output, option, P0)
      call write_Kijk  (output, option, boundary, Kijk)
      call write_Mij   (output, option, boundary, Mij)

      if (option%extrapolate) then

        write(iw,*)
        write(iw,'("Analyze> Start propagation")')

        ! Setup Boundary conditions 
        !
        call set_reflection(output, option, boundary, Kijk, Mij)
        call set_product   (output, option, boundary, Kijk, Mij)
        
        ! Extend timescale 
        !
        call reacdyn_tcf(output, option, boundary, Rij, P0, Kijk, Mij)

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
            write(iw,'("State determination failed. Please check state definition.")')
            stop
          end if

        end do
      end do 

    end subroutine get_state
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_f_init_id(option, nfile, init_id)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: nfile 
      integer,        intent(inout) :: init_id(option%nmol)

      ! I/O
      !
      integer :: io

      ! Dummy
      !
      integer :: ifile, imol


      ! Allocate
      !
      init_id = 0

      ! Read
      !
      call open_file(option%f_init_id, io, stat = 'old')
      do ifile = 1, nfile
        read(io,*) (init_id(imol), imol = 1, option%nmol)
      end do
      close(io)
      
!
    end subroutine read_f_init_id 
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

end module mod_analyze
!=======================================================================
