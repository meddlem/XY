program main
  use constants
  use initialize
  use markov
  use plotroutines
  use io
  implicit none
  
  ! variables:
  ! BE: beta*energy,
  ! dE: change in energy between new and initial config
  ! h: external field
  ! S: array containing Spins indexed as row, column

  real(dp), allocatable :: S(:,:,:), BE(:) 
  real(dp)              :: BK
  integer, allocatable  :: t(:)
  integer               :: runtime
  
  allocate(S(2,L,L),t(n_meas),BE(n_meas))
  
  call user_in(BK)
  call init_random_seed()
  call init_lattice(S)
  call animate_lattice('')
  
  call run_sim(S,BE,BK,t,runtime)
  
  call close_lattice_plot()
  call results_out(BK,BE(n_meas),runtime)
  call line_plot(real(t,dp),BE,'t','energy','','',1)
!  call line_plot(real(t,dp),real(m,dp),'t','magnetization','','',2)
  
!  deallocate(S,m,t,BE)
end program
