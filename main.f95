program main
  use constants
  use initialize
  use markov
  use plotroutines
  use io
  implicit none
  
  ! variables:
  ! BK: beta*coupling 
  ! BE: beta*energy 
  ! S: array containing spins vectors indexed as row, column

  real(dp), allocatable :: S(:,:,:), BE(:) 
  real(dp)              :: BK, h_mod
  integer, allocatable  :: t(:)
  integer               :: runtime
  
  allocate(S(2,L,L),t(n_meas),BE(n_meas))
  
  call user_in(BK)
  call init_random_seed()
  call init_lattice(S)
  call animate_lattice('')
  
  call run_sim(S,BE,BK,t,h_mod,runtime)
  
  call close_lattice_plot()
  call results_out(BK,BE(n_meas),h_mod,runtime)
  call line_plot(real(t,dp),BE,'t','energy','','',1)
  
  deallocate(S,t,BE)
end program
