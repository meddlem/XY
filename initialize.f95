module initialize
  use constants
  implicit none
  private
  public :: init_random_seed, init_lattice, init_energy

contains
  subroutine init_energy(BE,S,BJ,h)
    real(dp), intent(out) :: BE
    integer, intent(in) :: S(:,:)
    real(dp), intent(in) :: h, BJ
    integer :: i, j

    if (size(S,1) < 2) return !check
    
    BE = 0._dp ! initialze energy 
    
    do i = 2,L+1
      do j = 2,L+1
        BE = BE - BJ*S(i,j)*(S(i-1,j) + S(i+1,j) + S(i,j-1) + S(i,j+1))
      enddo
    enddo

    BE = 0.5_dp*BE ! account for double counting of pairs
    BE = BE - h*sum(S) ! add external field
  end subroutine

  subroutine init_lattice(S)
    integer, intent(out) :: S(:,:)
    real(dp), allocatable :: u(:,:)
    integer, allocatable :: S_tmp(:,:)

    allocate(u(L,L),S_tmp(L,L))
    S = 0

    call random_number(u)
    ! assign initial spins at random, corresponds to T=∞ 
    S_tmp = -1
    where (u>0.5_dp) S_tmp = 1
    S(2:L+1,2:L+1) = S_tmp
    deallocate(u)
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
