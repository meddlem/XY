program main
  use constants
  use initialize
  use markov
  implicit none

  real(dp) :: BJ = 10._dp, BE, dE, h = 0._dp
  integer :: i, j, S(L,L)
 
  ! dE_vals = 2._dp*BJ*(/(i,i=-4,4)/) ! possible values of dE
  ! BF_vals = exp(-dE_vals) ! possible values of boltzmann factor

  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE,S,BJ,h)

  do i=1,L
    write(*,'(100I3)')(S(i,j),j=1,L)
  enddo
  print *, 'initial energy:', BE
  
  do i=1,10*N
    call gen_config(S,dE,BJ,h)
    BE = BE + dE
  enddo

  do i=1,L
    write(*,'(100I3)')(S(i,j),j=1,L)
  enddo

  print *, 'energy:', BE
end program
