program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  ! variables:
  ! J: coupling 
  ! E: energy 
  ! L:  lattice side
  ! S:  array containing spins vectors indexed as row, column

  real(dp), allocatable :: S(:,:,:), t(:), E 
  real(dp)              :: J, h_mod, Chi
  integer               :: runtime, L
  
  call user_in(J,L)
  allocate(S(2,L,L),t(n_meas),E(n_meas))
  call init_random_seed()
  call init_lattice(L,S)
  
  call run_sim(S,E,J,t,h_mod,Chi,runtime)
  
  call results_out(J,t,E,h_mod,Chi,runtime)
  deallocate(S,t,E)
end program
