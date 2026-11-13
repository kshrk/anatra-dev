!=======================================================================
module mod_analyze
!=======================================================================
  use mod_util
  use mod_const
  use mod_grid3d
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_cv

  ! constants
  !
  real(8), parameter, private :: EPS = 1.0d-20

  ! subroutines
  !
  public  :: analyze
  private :: read_f_unperturbed_id
  private :: normalize_radial

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fhead_out, fout

      ! Local
      !
      integer :: nfile, nstate, ndim, ncol
      integer :: xyzcol(MaxDim), ng3(MaxDim), gind(MaxDim)
      real(8) :: origin(MaxDim), del(MaxDim)
      real(8) :: dv, val, kT, psum
      logical :: is_end, is_inside, is_detected

      type(s_cv)     :: cv
      type(s_func3d) :: gav, fel 

      ! Dummy
      !
      integer :: ifile, istep, idim, is, js, id
      integer :: ig

      ! Arrays
      !
      type(s_func3d), allocatable :: g(:)
      integer,        allocatable :: unperturbed_ids(:)
      integer,        allocatable :: nfile_each_state(:)
      real(8),        allocatable :: nsample(:)
      real(8),        allocatable :: pop(:, :), popav(:)


      ! Setup
      !
      nfile  = input%ncv
      nstate = option%nstate
      ndim   = option%ndim
      xyzcol = option%xyzcol
      ng3    = option%ng3
      del    = option%del
      origin = option%origin
      kT     = Boltz * option%temperature

      ! Check
      !
      if (ndim < 1 .or. ndim > MaxDim) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("ndim should be between 1 and ", i0)') MaxDim
        stop
      end if

      if (option%nmol > 1) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("nmol > 1 is currently not supported.")')
        stop
      end if

      if (.not. option%is_funp_specified) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("f_unperturbed_id should be specified.")')
        stop
      end if

      if (option%use_radial) then
        if (option%radcol > ndim) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("radcol should be less than or equal to ndim.")')
          stop
        end if
      end if

      ncol = maxval(xyzcol(1:ndim))

      ! Calculate measure dv = dx_1 * dx_2 * ...
      !
      dv = 1.0d0
      do idim = 1, ndim
        dv = dv * del(idim)
      end do

      ! Allocate
      !
      allocate(unperturbed_ids(nfile))
      allocate(nfile_each_state(nstate))
      allocate(nsample(nstate))
      allocate(g(nstate))
      allocate(pop(nstate, nstate))
      allocate(popav(nstate))

      nfile_each_state = 0
      nsample          = 0.0d0
      pop              = 0.0d0
      popav            = 0.0d0

      do is = 1, nstate
        call setup_func3d(ng3, del, origin, g(is))
      end do
      call setup_func3d(ng3, del, origin, gav)
      call setup_func3d(ng3, del, origin, fel)

      ! Read unperturbed-state ID for each file
      !
      write(iw,*)
      write(iw,'("Analyze> Read f_unperturbed_id file")')
      call read_f_unperturbed_id(option, nfile, unperturbed_ids)

      do ifile = 1, nfile
        id = unperturbed_ids(ifile)
        if (id < 1 .or. id > nstate) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("unperturbed state ID of file ", i0, " (", i0, ") &
                     &is out of range [1, nstate].")') ifile, id
          stop
        end if
        nfile_each_state(id) = nfile_each_state(id) + 1
      end do

      do is = 1, nstate
        write(iw,'("State ", i5, " : # of files = ", i0)') is, nfile_each_state(is)
      end do
      write(iw,'(">> Done")')

      ! Generate histogram for each unperturbed state
      !
      write(iw,*)
      write(iw,'("Analyze> Calculate Histogram")')

      do ifile = 1, nfile

        id = unperturbed_ids(ifile)

        write(iw,'("Read CV file: ", 2x, a, 2x, "(unperturbed state ", i0, ")")') &
          trim(input%fcv(ifile)), id

        is_end = .false.
        call open_file(input%fcv(ifile), io, stat = 'old')
        do while (.not. is_end)

          if (allocated(cv%x)) then
            call deallocate_cv(cv)
          end if

          call read_cv_split(io, ncol, '#SPLIT', cv, is_end)

          if (cv%nstep == 0) cycle

          do istep = 1, cv%nstep

            is_inside = .true.
            do idim = 1, ndim
              val        = cv%data(xyzcol(idim), istep)
              !if (val < option%state_def(1, idim, id) - del(idim) .or. &
              !    val > option%state_def(2, idim, id) + del(idim)) then
              !  is_inside = .false.
              !end if

              gind(idim) = nint((val - origin(idim)) / del(idim)) + 1

              if (gind(idim) < 1 .or. gind(idim) > ng3(idim)) then
                is_inside = .false.
              end if
            end do

            do is = 1, nstate

              is_detected = .true.
              do idim = 1, ndim
                val = cv%data(xyzcol(idim), istep)
                if (val <  option%state_def(1, idim, is) .or. &
                    val >= option%state_def(2, idim, is)) then
                  is_detected = .false.
                end if
              end do

              if (is_detected) then
                pop(is, id) = pop(is, id) + 1.0d0
              end if

            end do

            if (ndim < MaxDim) gind(ndim + 1:MaxDim) = 1

            if (is_inside) then
              g(id)%data(gind(1), gind(2), gind(3)) &
                = g(id)%data(gind(1), gind(2), gind(3)) + 1.0d0
              nsample(id) = nsample(id) + 1.0d0
            end if

          end do

        end do
        close(io)

      end do

      ! Normalize: int g(x) dx = 1 for each state
      !
      popav = 0.0d0
      do is = 1, nstate
        if (nsample(is) > 0.0d0) then
          g(is)%data = g(is)%data / (nsample(is) * dv)
        else
          write(iw,'("Analyze> Warning.")')
          write(iw,'("No sample was found inside the grid for state ", i0)') is
        end if

        psum       = sum(pop(:, is))
        if (psum > 0.0d0) then
          pop(:, is) = pop(:, is) / psum 
        end if

        gav%data = gav%data + option%state_weight(is) * g(is)%data
        popav(:) = popav(:) + option%state_weight(is) * pop(:, is)

        !! Normalization for radial coordinate
        !!
        !if (option%use_radial) then
        !  call normalize_radial(option, g(is))
        !end if
      end do

      psum  = sum(popav)
      popav = popav / psum 

      ! Normalization for radial coordinate
      !
      if (option%use_radial) then
        call normalize_radial(option, gav)
      end if

      write(iw,'(">> Finished")')

      ! Output
      !
      write(iw,*)
      write(iw,'("Analyze> Write histograms")')
      !do is = 1, nstate
      !  write(fhead_out,'(a,"_g_state",i3.3)') trim(output%fhead), is
      ! 
      !   if (option%ndim == 1) then
      !     call generate_1dmap_gnuplot(fhead_out, g(is))
      !   else if (option%ndim == 2) then
      !     call generate_2dmap_matplotlib(fhead_out, g(is))
      !     call generate_2dmap_gnuplot(fhead_out, g(is))
      !   else if (option%ndim == 3) then
      !     call write_dx(fhead_out, g(is))
      !   end if
      !end do

      write(fhead_out,'(a,"_g")') trim(output%fhead)
      !if (option%ndim == 1) then
      !  call generate_1dmap_gnuplot(fhead_out, gav)
      !else if (option%ndim == 2) then
      !  call generate_2dmap_matplotlib(fhead_out, gav)
      !  call generate_2dmap_gnuplot(fhead_out, gav)
      !else if (option%ndim == 3) then
      !  call write_dx(fhead_out, gav)
      !end if

      write(fhead_out,'(a,"_pmf")') trim(output%fhead)
      call generate_fel(kT, gav, fel)
      if (option%ndim == 1) then
        call generate_1dmap_gnuplot(fhead_out, fel)
      else if (option%ndim == 2) then
        call generate_2dmap_matplotlib(fhead_out, fel)
        call generate_2dmap_gnuplot(fhead_out, fel)
      else if (option%ndim == 3) then
        call write_dx(fhead_out, fel)
      end if

      write(fout, '(a,".pop")') trim(output%fhead) 
      call open_file(fout, io)
      do is = 1, nstate
        write(io,'(i10, 2x, e15.7)') is, popav(is)
      end do
      close(io)
      write(iw,'(">> Finished")')

      call deallocate_cv(cv)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_f_unperturbed_id(option, nfile, unperturbed_ids)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: nfile
      integer,        intent(inout) :: unperturbed_ids(nfile)

      ! I/O
      !
      integer :: io

      ! Dummy
      !
      integer :: ifile


      ! Note: only the first column of each line is used
      !       (compatible with the IEPDYN format: "unperturbed_id  use_for_Rij")
      !
      unperturbed_ids = 0

      call open_file(option%f_unperturbed_id, io, stat = 'old')
      do ifile = 1, nfile
        read(io,*) unperturbed_ids(ifile)
      end do
      close(io)

    end subroutine read_f_unperturbed_id
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_radial(option, g)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_func3d), intent(inout) :: g

      integer :: igr, rc
      real(8) :: dr, r, fact, fourpi


      rc     = option%radcol
      dr     = g%del(rc)
      fourpi = 4.0d0 * PI

      do igr = 1, g%ng3(rc)
        r    = g%origin(rc) + (igr - 1) * dr
        fact = fourpi * r * r
        if (r < EPS) then
          fact = fourpi / 3.0d0 * dr * dr
        end if

        if (rc == 1) then
          g%data(igr, :, :) = g%data(igr, :, :) / fact
        else if (rc == 2) then
          g%data(:, igr, :) = g%data(:, igr, :) / fact
        else if (rc == 3) then
          g%data(:, :, igr) = g%data(:, :, igr) / fact
        end if
      end do

    end subroutine normalize_radial
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_fel(kT, g, fel)
!-----------------------------------------------------------------------
      implicit none

      real(8),        intent(in)    :: kT
      type(s_func3d), intent(in)    :: g
      type(s_func3d), intent(inout) :: fel

      integer :: igx, igy, igz
      real(8) :: val, fmin


      !gmax    = maxval(g%data)
      do igz = 1, g%ng3(3)
        do igy = 1, g%ng3(2)
          do igx = 1, g%ng3(1)
            !val = g%data(igx, igy, igz) / gmax
            val = g%data(igx, igy, igz)
            if (val <= EPS) then
              fel%data(igx, igy, igz) = -kT*log(EPS) 
            else
              fel%data(igx, igy, igz) = -kT*log(val)
            end if
          end do
        end do
      end do

      fmin     = minval(fel%data)
      fel%data = fel%data - fmin 

    end subroutine generate_fel
!-----------------------------------------------------------------------


end module mod_analyze
!=======================================================================
