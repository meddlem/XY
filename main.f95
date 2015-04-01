program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  ! variables:
  ! BK: beta*coupling 
  ! BE: beta*energy 
  ! L:  lattice side
  ! S:  array containing spins vectors indexed as row, column

  real(dp), allocatable :: S(:,:,:), t(:), BE(:) 
  real(dp)              :: BK, h_mod, Xi
  integer               :: runtime, L
  
  call user_in(BK,L)
  allocate(S(2,L,L),t(n_meas),BE(n_meas))
  call init_random_seed()
  call init_lattice(L,S)
  
  call run_sim(S,BE,BK,t,h_mod,Xi,runtime)
  
  call results_out(BK,t,BE,h_mod,Xi,runtime)
  deallocate(S,t,BE)
end program
