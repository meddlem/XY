module initialize
  use constants
  implicit none
  private
  public :: init_lattice, init_energy

contains
  subroutine init_energy(BE,S,BJ)
    real(dp), intent(out) :: BE
    integer, intent(in) :: S(:,:)
    real(dp), intent(in) :: BJ
    integer :: i, j

    BE = 0._dp ! initialze energy 
    if (size(S,1)<2) return
    
    ! internal spins
    do i = 2,L-1
      do j = 2,L-1
        BE = BE + BJ*S(i,j)*(S(i-1,j) + S(i+1,j) + S(i,j-1) + S(i,j+1))
      enddo
    enddo

    ! top boundary
    do i = 2,L-1
      BE = BE + BJ*S(1,i)*(S(1,i-1)+S(1,i+1)+S(2,i))
    enddo
    ! bottom boundary
    do i = 2,L-1
      BE = BE + BJ*S(L,i)*(S(L,i-1)+S(L,i+1)+S(L-1,i))
    enddo
    ! left boundary
    do j = 2,L-1
      BE = BE + BJ*S(j,1)*(S(j-1,1)+S(j+1,1)+S(j,2))
    enddo
    ! right boundary
    do j = 2,L-1
      BE = BE + BJ*S(j,L)*(S(j-1,L)+S(j+1,L)+S(j,L-1))
    enddo

    ! add interactions of corner spins
    BE = BE + BJ*S(1,1)*(S(1,2)+S(2,1)) 
    BE = BE + BJ*S(1,L)*(S(1,L-1)+S(2,L)) 
    BE = BE + BJ*S(L,1)*(S(L,2)+S(L-1,1)) 
    BE = BE + BJ*S(L,L)*(S(L-1,L)+S(L,L-1))

    BE = 0.5_dp*BE ! account for double counting of pairs
  end subroutine

  subroutine init_lattice(S)
    integer, intent(out) :: S(L,L)
    integer :: i, j
    real(dp) :: u(L,L)

    call init_random_seed()
    call random_number(u)
    ! assign spins based on uniform random number u 
    do i = 1,L
      do j = 1,L
        if (u(i,j)>0.5_dp) then
          S(i,j) = 1
        else
          S(i,j) = -1
        endif
      enddo
    enddo
  end subroutine 

  ! initialize random seed, taken from ICCP github
  subroutine init_random_seed()
    integer, allocatable :: seed(:)
    integer :: i, m, un, istat, dtime(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = m)
    allocate(seed(m))
    open(newunit=un, file="/dev/urandom", access="stream",&
      form="unformatted", action="read", status="old", &
      iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      call system_clock(count)
      if (count /= 0) then
        t = transfer(count, t)
      else
        call date_and_time(values=dtime)
        tms = (dtime(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
          + dtime(2) * 31_8 * 24 * 60 * 60 * 1000 &
          + dtime(3) * 24 * 60 * 60 * 60 * 1000 &
          + dtime(5) * 60 * 60 * 1000 &
          + dtime(6) * 60 * 1000 + dtime(7) * 1000 &
          + dtime(8)
        t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (m >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (m > 3) then
          seed(4:) = s + 37 * (/ (i, i = 0, m - 4) /)
        end if
      else
        seed = s + 37 * (/ (i, i = 0, m - 1 ) /)
      end if
    end if
    call random_seed(put=seed)
  end subroutine
end module 
