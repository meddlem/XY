program main
  use constants
  use initialize
  implicit none
  real(dp) :: BJ =  1._dp, BE
  integer :: i, j, S(L,L)
  
  call init_lattice(S)
  do i=1,L
    write(*,'(100I3)')(S(i,j),j=1,L)
  enddo
  call init_energy(BE,S,BJ)
  print *, BE
end program
