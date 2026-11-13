!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)       :: input
  type(s_output)      :: output
  type(s_option)      :: option 

  call show_title
  call show_usage
  call read_ctrl  (input, output, option)
  call analyze    (input, output, option)
  call termination("IEPDYN-FE Analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("========================================================")')
  write(6,*)
  write(6,'("                I  E  P  D  Y  N - F E")')
  write(6,*)
  write(6,'("   Integral-Equation formalism of Population DYNamics")')
  write(6,*)
  write(6,'("========================================================")')
  write(6,'("[Developer]")')
  write(6,'("Kento Kasahara (The Univ. of Osaka)")')
  write(6,*)
  write(6,'("[Reference]")')
  write(6,'("K. Kasahara, R. Okabe, C. A. Chang, T. Mori, and N. Matubayasi,")')
  write(6,'("IEPDYN: Integral-equation formalism of population dynamics,")')
  write(6,'("J. Chem. Phys., 164, 124112 (2026).")')
  write(6,*)

end subroutine show_title
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine show_usage
!-----------------------------------------------------------------------
  use mod_const

  implicit none

  character(len=MaxChar) :: f_ctrl


  call getarg(1, f_ctrl)

  if (trim(f_ctrl) == "-h") then
    write(iw,'("&input_param")')
    write(iw,'(" fcv      = ""cvdata1"" ""cvdata2"" ! Time-series CV data")')
    write(iw,'(" flist_cv = ""list"" ! CV-file list (Either fcv or flist_cv should be specified)")')
    write(iw,'("/")')
    write(iw,'("&output_param")')
    write(iw,'(" fhead         = ""output""  ! header of output file")')
    write(iw,'("/")')
    write(iw,'("&option_param")')
    write(iw,'(" nstate      =  4             ! # of states")')
    write(iw,'(" ndim        =  1             ! # of dimensions (up to 3)")')
    write(iw,'(" nmol        =  1             ! # of target molecules (typically 1)")')
    write(iw,'(" temperature =  298.0d0       ! temperature [K]")')
    write(iw,*)
    write(iw,'(" xyzcol      =  1             ! column ID corresponding to x, y, z components in cv file")')
    write(iw,'("                              ! (please ignore the first column in cv file)")')
    write(iw,'(" ng3         =  10            ! # of grids for each cv")')
    write(iw,'(" del         =  0.1d0         ! grid spacing for each cv")')
    write(iw,'(" origin      =  0.0d0         ! origin for each cv")')
    write(iw,*)
    write(iw,'(" use_radial  =  .false.       ! use radial coordinate or not")')
    write(iw,'(" radcol      =  0             ! column ID of radial coordinate")')
    write(iw,'("                              ! (should be specified if use_radial = .true.)")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&state")')
    write(iw,'("! (lower, upper) of each dimension, followed by weight")')
    write(iw,'("-500.0  -16.0  0.5  ! State 1 (weight of 0.5)")')
    write(iw,'(" -16.0  -15.0  0.5  ! State 2 (weight of 0.5)")')
    write(iw,'(" -15.0   -2.0  0.0  ! State 3")')
    write(iw,'("  -2.0  500.0  0.0  ! State 4")')
    write(iw,'("/")')
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
