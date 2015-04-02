program main
  use constants
  use initialize
  use markov
  use io
  implicit none
  
  ! variables:
  ! BJ: coupling 
  ! BE: energy 
  ! L:  lattice side
  ! S:  array containing spins vectors indexed as row, column

  real(dp), allocatable :: S(:,:,:), t(:), BE(:)
  real(dp)              :: BJ, h_mod, Chi
  integer               :: runtime, L
  
  call user_in(BJ,L)
  allocate(S(2,L,L),t(n_meas),BE(n_meas))
  call init_random_seed()
  call init_lattice(L,S)
  
  call run_sim(S,BE,BJ,t,h_mod,Chi,runtime)
  
  call results_out(BJ,t,BE,h_mod,Chi,runtime)
  deallocate(S,t,BE)
end program
