program main
  use constants
  use initialize
  use markov
  implicit none

  real(dp) :: BJ =  1._dp, BE, dE, h = 0._dp
  integer :: i, j, S(L,L)
  
  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE,S,BJ,h)

  do i=1,L
    write(*,'(100I3)')(S(i,j),j=1,L)
  enddo
  print *, 'initial energy:', BE
 
  call gen_config(S,dE,BJ,h)
  BE = BE + dE

  do i=1,L
    write(*,'(100I3)')(S(i,j),j=1,L)
  enddo

  print *, 'energy:', BE
end program
