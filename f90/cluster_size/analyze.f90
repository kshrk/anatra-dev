!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_ctrl
  use mod_xtcio
  use mod_dcdio
  use mod_netcdfio
  use mod_traj
  use mod_com
  use xdr, only: xtcfile

  ! subroutines
  !
  public  :: analyze
  private :: union
  private :: find

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),    intent(in)    :: input 
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_traj),     intent(inout) :: traj

      ! I/O
      !
      integer                :: io

      ! Local
      ! 
      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc
      type(s_netcdf)         :: nc 
      type(s_com)            :: com

      integer                :: nmol, natm
      integer                :: trajtype
      logical                :: is_end
      real(8)                :: dvec(3), dsq, rcsq, acc
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer                :: itraj, imol, jmol, isize, iroot
      integer                :: istep, istep_tot

      ! Arrays
      !
      integer, allocatable   :: parent(:)
      integer, allocatable   :: cluster_size(:)
      real(8), allocatable   :: distr(:)
      !real(8), allocatable   ::


      ! Setup
      !
      rcsq = option%rcut * option%rcut

      call get_com(option%mode,          &
                   traj,                 &
                   com,                  &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)
      nmol = com%nmol

      ! Allocate
      !
      allocate(parent(nmol),        &
               cluster_size(nmol),  &
               distr(nmol))

      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Start 
      !
      write(iw,*)
      write(iw,'("Analyze> Start")')
      istep_tot = 0
      distr     = 0.0d0
      do itraj = 1, input%ntraj

        call open_trajfile(input%ftraj(itraj), trajtype, io, dcd, xtc, nc)
        call init_trajfile(trajtype, io, dcd, xtc, nc, natm)

        is_end = .false.
        istep  = 0
        do while (.not. is_end)

          cluster_size = 0
          istep        = istep     + 1
          istep_tot    = istep_tot + 1

          call read_trajfile_oneframe(trajtype, io, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ", i0)') istep_tot
          end if

          call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj)

          call get_com(option%mode,         &
                       traj,                &
                       com,                 &
                       setup = .false.,     &
                       calc_coord = .true., & 
                       myrank = 1)


          ! Initialize 
          !
          do imol = 1, nmol
            parent(imol) = imol
          end do

          ! Union
          !
          do imol = 1, nmol - 1
            do jmol = imol + 1, nmol
              dvec(1:3) = com%coord(1:3, imol, 1) - com%coord(1:3, jmol, 1)
              dvec(1:3) = dvec(1:3) - anint(dvec(1:3)/traj%box(1:3, 1)) * traj%box(1:3, 1)
              dsq       = dot_product(dvec(1:3), dvec(1:3))

              if (dsq <= rcsq) then
                call union(imol, jmol, parent)
              end if
            end do
          end do

          ! Find 
          !
          do imol = 1, nmol
            iroot               = find(imol, parent)
            cluster_size(iroot) = cluster_size(iroot) + 1
          end do

          ! Calc. histogram 
          !
          do iroot = 1, nmol
            isize        = cluster_size(iroot)

            if (isize == 0) cycle
            distr(isize) = distr(isize) + 1.0d0 
          end do

        end do

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do

      ! Normalize 
      !
      !distr(:) = distr(:) / (dble(istep_tot * nmol))
      distr(:) = distr(:) / (dble(istep_tot))
      distr(:) = distr(:) / sum(distr(:))

      ! Output 
      !
      write(fname, '(a, ".distr")') trim(output%fhead)
      call open_file(fname, io)
      acc = 0.0d0
      do isize = 1, nmol
        !acc = acc + distr(isize) * isize
        acc = acc + distr(isize)
        write(io, '(i0,2x,e15.7,2x,e15.7)') isize, distr(isize), acc
      end do 
      close(io) 


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    recursive function find(x, parent) result(ret)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)    :: x
      integer, intent(inout) :: parent(:)

      integer :: ret


      if (parent(x) == x) then
        ret = x
        return
      else
        parent(x) = find(parent(x), parent)
        ret       = parent(x)
      end if

    end function find
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine union(x, y, parent)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)    :: x
      integer, intent(in)    :: y
      integer, intent(inout) :: parent(:)

      ! Local
      !
      integer :: root_x, root_y


      root_x = find(x, parent)
      root_y = find(y, parent)

      if (root_x /= root_y) then
        parent(root_x) = root_y
      end if

    end subroutine union
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
