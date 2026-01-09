!-----------------------------------------------------------------------
    subroutine reacdyn_pint(output, option, boundary, Rij, P0, Kijk, Mij)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      real(8),          intent(in) :: Rij(0:option%nt_range, &
                                          option%nstate,     &
                                          option%nstate)
      real(8),          intent(in) :: P0(0:option%nt_range,  &
                                         option%nstate)
      real(8),          intent(in) :: Kijk(0:option%nt_range, &
                                          option%nstate,      &
                                          -boundary%nboundary:boundary%nboundary)
      real(8),          intent(in) :: Mij(0:option%nt_range,  &
                                         -boundary%nboundary:boundary%nboundary)

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range, nboundary, nbt
      real(8)                :: dt
      character(len=MaxChar) :: fname

      integer                :: thread_id

      ! Dummy
      !
      integer :: is, js, ks, is1, is2, js1, js2
      integer :: iu, ju, ku
      integer :: ib, jb, inflx, jnflx, id
      integer :: lda, ldb, info
      real(8) :: psum

      ! Arrays
      !
      real(8), allocatable :: Qint(:)
      real(8), allocatable :: Kuu(:, :, :), Mu(:, :), Ru(:, :)
      real(8), allocatable :: Xuu(:, :), Pint(:)
      integer, allocatable :: map(:, :), mapb(:), mapb_inv(:)
      integer, allocatable :: ipiv(:)


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      dt        = option%dt_out
      nboundary = boundary%nboundary
      nbt       = nboundary * 2

      allocate(Mu(0:nt_range, nbt))
      allocate(Kuu(0:nt_range, nbt, nbt))
      allocate(Ru(0:nt_range, nbt))
      allocate(Xuu(nbt, nbt))
      allocate(Pint(nstate))
      allocate(Qint(nbt))
      allocate(ipiv(nbt)) 
      allocate(map(2, nbt))
      allocate(mapb(nbt))
      allocate(mapb_inv(-nboundary:nboundary))
      Mu       = 0.0d0
      Kuu      = 0.0d0
      Ru       = 0.0d0
      Xuu      = 0.0d0
      Pint     = 0.0d0
      Qint     = 0.0d0
      map      = 0
      mapb     = 0
      mapb_inv = 0

      iu  = 0
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        iu  = iu + 1
        is1 = boundary%b2P(1, ib)
        is2 = boundary%b2p(2, ib)
        map(1, iu)          = is1
        map(2, iu)          = is2
        mapb(iu)            = ib
        mapb_inv(ib)        = iu

        Mu(0:nt_range, iu) = Mij(0:nt_range, ib)
        Ru(0:nt_range, iu) = Rij(0:nt_range, is2, is1)

        ju = 0
        do jb = -nboundary, nboundary
          if (jb == 0)  cycle
          ju  = ju + 1
          js1 = boundary%b2p(1, jb)
          js2 = boundary%b2p(2, jb)

          if (is1 == js2) then
            Kuu(0:nt_range, iu, ju) = Kijk(0:nt_range, is2, jb) 
          end if 
        end do
      end do

      ! Reflecting boundary
      !
      if (option%use_reflection_state) then
        do iu = 1, nbt
          ib = mapb(iu)
          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)
          if (boundary%conv_direc(ib)) then
            jb = -ib
            ju = mapb_inv(jb)
            Ru(0:nt_range, ju) = Rij(0:nt_range, is2, is1)
            Ru(0:nt_range, iu) = 0.0d0

            Kuu(0:nt_range, ju, :) = 0.0d0
            do ku = 1, nbt
              Kuu(0:nt_range, ju, ku) = Kuu(0:nt_range, iu, ku)
              Kuu(0:nt_range, iu, ku) = 0.0d0
            end do 
          end if
        end do
      end if

      do ju = 1, nbt
        Xuu(ju, ju) = 1.0d0
        do iu = 1, nbt
          Xuu(iu, ju) = Xuu(iu, ju) - sum(Kuu(0:nt_range, iu, ju)) * dt
        end do
      end do

      do iu = 1, nbt
        Qint(iu) = sum(Ru(0:nt_range, iu)) * dt
      end do
      lda  = nbt
      ldb  = nbt
      ipiv = 0
      call dgesv(nbt, 1, Xuu, lda, ipiv, Qint, ldb, info) 

      do is = 1, nstate 
        if (option%is_dissoc(is)) then
          Pint(is) = 0.0d0
          cycle       
        end if

        Pint(is) = sum(P0(0:nt_range, is)) * dt
        do iu = 1, nbt
          if (map(2, iu) /= is) cycle
          Pint(is) = Pint(is) + sum(Mu(0:nt_range, iu)) * dt * Qint(iu)  
        end do
      end do

      ! Integrated values 
      !
      write(fname,'(a,".int")') trim(output%fhead)
      call open_file(fname, io)
      psum = 0.0
      do is = 1, nstate
        if (option%is_initial(is)) then
          psum = psum + Pint(is)
        end if
      end do

      write(io,'("Pint Total  ", e15.7)') psum
      do is = 1, nstate
        if (option%is_dissoc(is)) cycle
        write(io,'("Pint ", i5, 2x, e15.7)') is, Pint(is)
      end do

      close(io)

      ! Deallocate memory
      !
      
    end subroutine reacdyn_pint
!-----------------------------------------------------------------------
