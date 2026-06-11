!-------------------------------------------------------------------------------
!
!  Module   mod_fft
!
!    Module for using MKL-FFT library   
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

include 'mkl_dfti.f90'

module mod_fft
  use MKL_DFTI 

  ! constants
  !
  integer, parameter :: FFT_SignForward  = 1
  integer, parameter :: FFT_SignBackward = 2 

  ! structures
  !
  type s_fftinfo
    type(Dfti_Descriptor), pointer :: desc_forward
    type(Dfti_Descriptor), pointer :: desc_backward
    integer    :: ng3(3)
  end type s_fftinfo

  ! variables
  !
  integer :: fftsize(3)

  ! subroutines
  !
  public :: fft_init
  public :: fft_r2c
  public :: fft_c2r

  contains
!-----------------------------------------------------------------------
      subroutine fft_init(fftinfo, funck, funcm)
!-----------------------------------------------------------------------
      implicit none
!
      type(s_fftinfo),     intent(inout) :: fftinfo
      real(8),             intent(inout) :: funck(fftinfo%ng3(1),       &
                                                  fftinfo%ng3(2),       &
                                                  fftinfo%ng3(3))
      complex(kind(0d0)),  intent(inout) :: funcm(fftinfo%ng3(1)/2 + 1, &
                                                  fftinfo%ng3(2),       &
                                                  fftinfo%ng3(3))

      integer :: igx, igy, igz, igr, igk
      integer :: ngx, ngy, ngz, ngr, ngk
      integer :: statf, statb
      integer :: strides(4)

      
      ! setup module variables
      !
      fftsize(1:3) = fftinfo%ng3(1:3)

      
      ngx = fftinfo%ng3(1)
      ngy = fftinfo%ng3(2)
      ngz = fftinfo%ng3(3)

      strides(1) = 0
      strides(2) = 1
      strides(3) = ngx / 2 + 1
      strides(4) = strides(3) * ngy

      ! Forward
      !
      statf = DftiCreateDescriptor(fftinfo%desc_forward,  &
                                   DFTI_DOUBLE,           &
                                   DFTI_REAL,             &
                                   3,                     &
                                   fftinfo%ng3)

      statf = DftiSetValue(fftinfo%desc_forward,          &
                           DFTI_PLACEMENT,                &
                           DFTI_NOT_INPLACE)
   
      statf = DftiSetValue(fftinfo%desc_forward,          &
                           DFTI_CONJUGATE_EVEN_STORAGE,   &
                           DFTI_COMPLEX_COMPLEX)

      statf = DftiSetValue(fftinfo%desc_forward,          &
                           DFTI_OUTPUT_STRIDES,           &
                           strides)

      statf = DftiSetValue(fftinfo%desc_forward,          &
                          DFTI_PACKED_FORMAT,             &
                          DFTI_CCE_FORMAT)

      statf = DftiCommitDescriptor(fftinfo%desc_forward)

      ! Backward
      !
      statb = DftiCreateDescriptor(fftinfo%desc_backward, &
                                   DFTI_DOUBLE,           &
                                   DFTI_REAL,             &
                                   3,                     &
                                   fftinfo%ng3)

      statb = DftiSetValue(fftinfo%desc_backward,         &
                           DFTI_PLACEMENT,                &
                           DFTI_NOT_INPLACE)
   
      statb  = DftiSetValue(fftinfo%desc_backward,        &
                           DFTI_CONJUGATE_EVEN_STORAGE,   &
                           DFTI_COMPLEX_COMPLEX)

      statb = DftiSetValue(fftinfo%desc_backward,         &
                           DFTI_INPUT_STRIDES,            &
                           strides)

      statb = DftiSetValue(fftinfo%desc_backward,         &
                          DFTI_PACKED_FORMAT,             &
                          DFTI_CCE_FORMAT)

      statb = DftiCommitDescriptor(fftinfo%desc_backward)

      end subroutine fft_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fft_r2c(fftinfo, in, out)
!-----------------------------------------------------------------------
      use MKL_DFTI 
      implicit none 

      type(s_fftinfo),     intent(in)  :: fftinfo
      real(8),             intent(in)  :: in(*)
      complex(kind(0d0)),  intent(out) :: out(*)

      integer :: stat

      stat = DftiComputeForward(fftinfo%desc_forward, in, out)

      end subroutine fft_r2c
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fft_c2r(fftinfo, in, out)
!-----------------------------------------------------------------------
      use MKL_DFTI 
      implicit none 

      type(s_fftinfo),     intent(in)  :: fftinfo
      complex(kind(0d0)),  intent(in)  :: in(*)
      real(8),             intent(out) :: out(*)

      integer :: stat

      stat = DftiComputeBackward(fftinfo%desc_backward, in, out) 
!
      end subroutine fft_c2r
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fft_cleanup(fftinfo)
!-----------------------------------------------------------------------
      implicit none 

      type(s_fftinfo),     intent(in)  :: fftinfo

      integer :: stat

      stat = DftiFreeDescriptor(fftinfo%desc_forward)
      stat = DftiFreeDescriptor(fftinfo%desc_backward)
!
      end subroutine fft_cleanup
!----------------------------------------------------------------------

end module mod_fft
