!-----------------------------------------------------------------------
    subroutine smooth_RK(option, b, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: b
      type(s_func),     intent(inout) :: f

      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, js, is1, is2, ib, jb, inflx
      real(8) :: state_sum, rsum, ksum

      ! Arrays
      !
      real(8), allocatable :: wrk(:) 


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      dt        = option%dt_out
      nboundary = b%nboundary

      if (.not. allocated(wrk)) then
        allocate(wrk(0:nt_range))
      end if

      ! Rij
      !
      do js = 1, nstate
        do is = 1, nstate
          rsum = sum(f%R(:, is, js))
          if (rsum < 1.0d-10) cycle
          call smooth_histogram(nt_range + 1,   &
                                f%R(0, is, js), &
                                wrk(0),         &
                                option%smooth_order)
          f%R(0:nt_range, is, js) = wrk(0:nt_range)
        end do
      end do

      ! Kijk
      !
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        is1 = b%b2p(1, ib)
        is2 = b%b2p(2, ib)
        do inflx = 1, b%n_influx_boundary(is1)
          jb    = b%influx_boundary(inflx, is1)
          js    = b%b2p(1, jb)
          ksum  = sum(f%K(:, is2, jb))
          if (ksum < 1.0d-10) cycle
          call smooth_histogram(nt_range + 1,    &
                                f%K(0, is2, jb), &
                                wrk(0),          &
                                option%smooth_order)
          f%K(0:nt_range, is2, jb) = wrk(0:nt_range)
        end do 
      end do 

    end subroutine smooth_RK 
!-----------------------------------------------------------------------

subroutine smooth_histogram(n, x, y, nk)

    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: nk          ! odd kernel size (3,5,7,...)

    real(8), intent(in)  :: x(n)
    real(8), intent(out) :: y(n)

    integer :: i, j
    integer :: idx
    integer :: m

    real(8), allocatable :: kernel(:)

    real(8) :: value
    real(8) :: sumw
    real(8) :: orig_int
    real(8) :: smooth_int
    real(8) :: scale

    !--------------------------------------------------
    ! Check kernel size
    !--------------------------------------------------

    if (mod(nk,2) == 0 .or. nk < 3) then
        write(*,*) "ERROR: nk must be an odd integer >= 3."
        stop
    endif

    m = nk/2

    allocate(kernel(nk))

    !--------------------------------------------------
    ! Generate Gaussian-approximation kernel
    ! using binomial coefficients
    !--------------------------------------------------

    kernel(1) = 1.d0

    do i = 2, nk
        kernel(i) = kernel(i-1) * dble(nk-i+1) / dble(i-1)
    enddo

    kernel = kernel / sum(kernel)

    !--------------------------------------------------
    ! Original rectangular integral
    !--------------------------------------------------

    orig_int = sum(x)

    !--------------------------------------------------
    ! Convolution
    !--------------------------------------------------

    do i = 1, n

        value = 0.d0
        sumw  = 0.d0

        do j = -m, m

            idx = i + j

            if (idx >= 1 .and. idx <= n) then
                value = value + kernel(j+m+1) * x(idx)
                sumw  = sumw  + kernel(j+m+1)
            endif

        enddo

        y(i) = value / sumw

    enddo

    !--------------------------------------------------
    ! Preserve rectangular integral
    !--------------------------------------------------

    smooth_int = sum(y)

    if (abs(smooth_int) > tiny(smooth_int)) then
        scale = orig_int / smooth_int
        y = y * scale
    endif

    deallocate(kernel)

end subroutine smooth_histogram

!!-----------------------------------------------------------------------
!    subroutine smooth_histogram(n, x, y)
!!-----------------------------------------------------------------------
!
!      implicit none
!
!      integer, intent(in) :: n
!
!      real(8), intent(in)  :: x(n)
!      real(8), intent(out) :: y(n)
!
!      integer :: i
!      integer :: j
!      integer :: idx
!
!      real(8), parameter :: kernel(5) = (/ &
!          1.d0/16.d0, &
!          4.d0/16.d0, &
!          6.d0/16.d0, &
!          4.d0/16.d0, &
!          1.d0/16.d0 /)
!
!      real(8) :: value
!      real(8) :: sumw
!      real(8) :: orig_int
!      real(8) :: smooth_int
!      real(8) :: scale
!
!      !-----------------------------------------
!      ! Original rectangular integral
!      !-----------------------------------------
!
!      orig_int = sum(x)
!
!      do i = 1, n
!
!          value = 0.0d0
!          sumw  = 0.0d0
!
!          do j = -2, 2
!
!              idx = i + j
!
!              if (idx >= 1 .and. idx <= n) then
!                  value = value + kernel(j+3) * x(idx)
!                  sumw  = sumw  + kernel(j+3)
!              endif
!
!          enddo
!
!          y(i) = value / sumw
!
!      enddo
!
!      !-----------------------------------------
!      ! Preserve rectangular integral
!      !-----------------------------------------
!
!      smooth_int = sum(y)
!
!      if (abs(smooth_int) > tiny(smooth_int)) then
!          scale = orig_int / smooth_int
!          y = y * scale
!      endif
!
!
!    end subroutine smooth_histogram
!!-----------------------------------------------------------------------
