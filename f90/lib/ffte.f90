!-------------------------------------------------------------------------------
!
!  Module   mod_fft
!
!    Module for using FFTE library by Daisuke Takahashi
!    This module was developed by Kento Kasahara
!
!    (c) Copyright 2026 Kento Kasahara. All rights reserved.
!
!-------------------------------------------------------------------------------

module mod_fft

  ! constants
  !
  integer, parameter :: FFT_Initialize   =  0
  integer, parameter :: FFT_SignForward  = -1
  integer, parameter :: FFT_SignBackward = +1

  ! structures
  !
  type s_fftinfo
    integer              :: ng3(3)
    real(8), allocatable :: a(:), b(:) 
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
      subroutine fft_init(fi, in, out)
!-----------------------------------------------------------------------
      implicit none
!
      type(s_fftinfo),    intent(inout) :: fi
      real(8),            intent(inout) :: in(fi%ng3(1), fi%ng3(2), fi%ng3(3))
      complex(kind(0d0)), intent(inout) :: out(fi%ng3(1)/2+1, fi%ng3(2), fi%ng3(3)) 

      integer :: igx, igy, igz, igr, igk
      integer :: ngx, ngy, ngz, ngr, ngk
     

      ngx = fi%ng3(1)
      ngy = fi%ng3(2)
      ngz = fi%ng3(3)

      if (.not. allocated(fi%a)) then
        allocate(fi%a((ngx+2)*ngy*ngz))
        allocate(fi%b((ngx+2)*ngy*ngz))
        fi%a = 0.0d0
        fi%b = 0.0d0
      end if

      call dzfft3d(fi%a, ngx, ngy, ngz, FFT_Initialize, fi%b)
      call zdfft3d(fi%a, ngx, ngy, ngz, FFT_Initialize, fi%b)

      end subroutine fft_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fft_r2c(fi, in, out)
!-----------------------------------------------------------------------
      implicit none 

      type(s_fftinfo),    intent(inout) :: fi
      real(8),            intent(inout) :: in(fi%ng3(1), fi%ng3(2), fi%ng3(3))
      complex(kind(0d0)), intent(inout) :: out(fi%ng3(1)/2+1, fi%ng3(2), fi%ng3(3)) 

      integer :: ngx, ngy, ngz
      integer :: igx, igy, igz, igxo, igxe, ig, igo, ige 
    
      ngx = fi%ng3(1)
      ngy = fi%ng3(2)
      ngz = fi%ng3(3)

      fi%a = 0.0d0
      fi%b = 0.0d0

      do igz = 1, ngz
        do igy = 1, ngy
          do igx = 1, ngx
            ig = igx + ngx * (igy - 1) + ngx * ngy * (igz - 1)
            fi%a(ig) = in(igx, igy, igz)
          end do
        end do
      end do

      call dzfft3d(fi%a, ngx, ngy, ngz, FFT_SignForward, fi%b)

      do igz = 1, ngz
        do igy = 1, ngy
          do igx = 1, ngx/2 + 1

            igxo = 2 * igx - 1
            igxe = 2 * igx

            igo = igxo + (ngx+2) * (igy - 1) + (ngx+2) * ngy * (igz - 1)
            ige = igxe + (ngx+2) * (igy - 1) + (ngx+2) * ngy * (igz - 1)

            out(igx, igy, igz) = dcmplx(fi%a(igo), fi%a(ige))
!            print*, igx, igy, igz, out(igx, igy, igz)

          end do
        end do
      end do

!      stop

      end subroutine fft_r2c
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fft_c2r(fi, in, out)
!-----------------------------------------------------------------------
      implicit none 

      type(s_fftinfo),    intent(inout) :: fi
      complex(kind(0d0)), intent(inout) :: in(fi%ng3(1)/2+1, fi%ng3(2), fi%ng3(3)) 
      real(8),            intent(inout) :: out(fi%ng3(1), fi%ng3(2), fi%ng3(3))

      integer :: ngx, ngy, ngz
      integer :: igx, igy, igz, igxo, igxe, ig, igo, ige 


      ngx = fi%ng3(1)
      ngy = fi%ng3(2)
      ngz = fi%ng3(3)

      fi%a = 0.0d0
      fi%b = 0.0d0
      do igz = 1, ngz
        do igy = 1, ngy
          do igx = 1, ngx/2 + 1

            igxo = 2*igx - 1
            igxe = 2*igx

            !igo = igxo + ngx * (igy - 1) + ngx * ngy * (igz - 1)
            !ige = igxe + ngx * (igy - 1) + ngx * ngy * (igz - 1)
            igo = igxo + (ngx+2) * (igy - 1) + (ngx+2) * ngy * (igz - 1)
            ige = igxe + (ngx+2) * (igy - 1) + (ngx+2) * ngy * (igz - 1)

            fi%a(igo) = dble (in(igx, igy, igz))
            fi%a(ige) = dimag(in(igx, igy, igz))

          end do
        end do
      end do

      call zdfft3d(fi%a, ngx, ngy, ngz, FFT_SignBackward, fi%b)

      do igz = 1, ngz
        do igy = 1, ngy
          do igx = 1, ngx
            ig = igx + ngx * (igy - 1) + ngx * ngy * (igz - 1)
            out(igx, igy, igz) = fi%a(ig) 
          end do
        end do
      end do

      out = out * dble(ngx*ngy*ngz)

    end subroutine fft_c2r
!----------------------------------------------------------------------

end module mod_fft
