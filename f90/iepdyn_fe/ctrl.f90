!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  implicit none

  ! constants
  !
  integer,      parameter, public :: MaxDim       =  3
  integer,      parameter, public :: MaxStates    =  100 
  integer,      parameter, public :: NotSpecified = -100

  ! structures
  !
  type :: s_option
    integer :: nmol                 = NotSpecified 
    integer :: ndim                 = NotSpecified 
    integer :: nstate               = NotSpecified

    ! Grid info.
    integer :: xyzcol(MaxDim)       = (/2, 3, 4/)
    integer :: ng3(MaxDim)          = 100
    real(8) :: del(MaxDim)          = (/0.1d0, 0.1d0, 0.1d0/) 
    real(8) :: origin(MaxDim)       = 0.0d0

    ! Radial coordinate
    logical :: use_radial           = .false.
    integer :: radcol               = 0

    ! File names 
    !
    character(len=MaxChar) :: f_unperturbed_id = ''

    ! for free-energy calculation 
    !
    real(8)              :: temperature  = 298.0d0

    ! prepared after reading namelists
    !
    logical              :: is_funp_specified = .false.

    real(8), allocatable :: state_def   (:, :, :)
    real(8), allocatable :: state_weight(:)
    real(8), allocatable :: state_weight_unnorm(:)

  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: read_ctrl_state
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(out) :: input
      type(s_output), intent(out) :: output
      type(s_option), intent(out) :: option 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: f_ctrl


      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      call open_file        (f_ctrl, io)
      call read_ctrl_input  (io, input)
      call show_input       (input)
      call read_ctrl_output (io, output)
      call show_output      (output)
      call read_ctrl_option (io, option)
      call read_ctrl_state  (io, option)
      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer, parameter :: ndim_max = 3 
!
      integer,          intent(in)  :: io
      type(s_option),   intent(out) :: option

      character(len=MaxChar) :: f_unperturbed_id = ''

      integer :: nmol                 = NotSpecified
      integer :: ndim                 = NotSpecified
      integer :: nstate               = NotSpecified
      integer :: xyzcol(MaxDim)       = (/1, 2, 3/)
      integer :: ng3(MaxDim)          = (/1, 1, 1/)
      real(8) :: del(MaxDim)          = (/1.0d0, 1.0d0, 1.0d0/) 
      real(8) :: origin(MaxDim)       = (/0.0d0, 0.0d0, 0.0d0/)
      logical :: use_radial           = .false.
      integer :: radcol               = 0
      real(8) :: temperature          = 300.0d0

      ! Dummy
      !
      integer :: i


      namelist /option_param/       &
        f_unperturbed_id,           &
        nmol,                       &
        ndim,                       &
        nstate,                     &
        xyzcol,                     &
        ng3,                        &
        del,                        &
        origin,                     &
        use_radial,                 &
        radcol,                     &
        temperature

      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')

      write(iw,'("f_unperturbed_id     = ", a)')      trim(f_unperturbed_id)
      write(iw,'("nmol                 = ", i0)')     nmol
      write(iw,'("ndim                 = ", i0)')     ndim
      write(iw,'("nstate               = ", i0)')     nstate
      write(iw,'("xyzcol               = ", 3(i0,2x))')        (xyzcol(i), i = 1, ndim)
      write(iw,'("ng3                  = ", 3(i0,2x))')        (ng3(i),    i = 1, ndim)
      write(iw,'("del                  = ", 3(f15.7,2x))')     (del(i),    i = 1, ndim)
      write(iw,'("origin               = ", 3(f15.7,2x))')     (origin(i), i = 1, ndim)
      write(iw,'("use_radial           = ", a)')      get_tof(use_radial)
      if (use_radial) then
        write(iw,'("radcol               = ", i0)')     radcol
      end if
      write(iw,'("temperature          = ", f20.10)') temperature

      ! Check
      !
      if (nmol == NotSpecified) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("nmol shouldb be specified.")')
        stop
      end if

      if (ndim == NotSpecified) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("ndim shouldb be specified.")')
        stop
      end if

      if (nstate == NotSpecified) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("nstate should be specified.")')
        stop
      end if

      if (use_radial) then
        if (radcol == 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("radcol (column id of radial coordinate) &
                     &should be specified")')
          stop
        end if
      end if

      option%is_funp_specified = .false.
      if (trim(f_unperturbed_id) /= '') then
        option%is_funp_specified = .true.
      end if

      option%f_unperturbed_id     = f_unperturbed_id
      option%nmol                 = nmol
      option%ndim                 = ndim
      option%nstate               = nstate
      option%temperature          = temperature

      !   Grid
      option%xyzcol               = 0       ! Initialize
      option%ng3                  = 1       ! Initialize
      option%del                  = 1.0d0   ! Initialize

      option%xyzcol(1:ndim)       = xyzcol(1:ndim)
      option%ng3(1:ndim)          = ng3(1:ndim)
      option%del(1:ndim)          = del(1:ndim)
      option%origin               = origin

      !   Radial
      option%use_radial           = use_radial
      option%radcol               = radcol

      ! Memory allocation
      !
      allocate(option%state_def(2, option%ndim, option%nstate))
      allocate(option%state_weight(option%nstate))
      allocate(option%state_weight_unnorm(option%nstate))


    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_state(io, option)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: io
      type(s_option), intent(inout) :: option

      ! Local
      !
      character(len=MaxChar) :: line
      integer                :: nstate
      real(8)                :: weight_sum

      ! Dummy
      !
      integer :: istate, idm, imm, ierr


      ! Setup
      !
      nstate   = option%nstate

      call seek_line(io, '&state', ierr)

      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_State> Error.")')
        write(iw,'("&state section is not found.")')
        stop
      end if

      istate                     = 0
      option%state_def           = 0.0d0
      option%state_weight        = 0.0d0
      option%state_weight_unnorm = 0.0d0
      do while (.true.)
        read(io,'(a)', end = 100) line

        if (line(1:1) /= "#" .and. line(1:1) /= "!") then
          istate = istate + 1

          read(line,*) ((option%state_def(imm, idm, istate), imm = 1, 2), &
                         idm = 1, option%ndim), option%state_weight(istate)
        end if

        if (istate == option%nstate) exit

      end do

100   continue

      write(iw,*)
      write(iw,'("Read_Ctrl_State>  States are defined as")')
      do istate = 1, option%nstate
        write(iw,'("State ", i5)') istate

        do idm = 1, option%ndim
          write(iw,'("  component ", i5, " : ", f15.7, " <===> ", f15.7)') &
            idm,                              &
            option%state_def(1, idm, istate), &
            option%state_def(2, idm, istate)
        end do

        write(iw,*)

      end do

      ! Weight
      ! (all states are treated as initial states)
      !
      weight_sum = sum(option%state_weight(1:nstate))
      option%state_weight_unnorm = option%state_weight
      option%state_weight        = option%state_weight / weight_sum 

      write(iw,*)
      write(iw,'("Read_Ctrl_State> Weight for each state")')
      do istate = 1, nstate
        write(iw,'(i5,f15.7)') istate, option%state_weight(istate)
      end do
      write(iw,*)

    end subroutine read_ctrl_state
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_input(input)
!-----------------------------------------------------------------------
      implicit none

      type(s_input), intent(in) :: input

      ! Dummy
      !
      integer :: i

      ! Check
      !
      if (input%ncv == 0) then
        write(iw,'("Error. fcv should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Input section parameters")')
      do i = 1, input%ncv
        write(iw,'("fcv", 3x, i0, 3x, " = ", a)') i, trim(input%fcv(i)) 
      end do


    end subroutine show_input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_output(output)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output 


      ! Check
      !
      if (trim(output%fhead) == '') then
        write(iw,'("Error. fhead should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Output section parameters")')
      write(iw,'("fhead          = ", a)') trim(output%fhead)


    end subroutine show_output
!-----------------------------------------------------------------------

end module
!=======================================================================
