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

  real(dp), allocatable :: BE(:)
  real(dp)              :: BJ, h
  integer, allocatable  :: S(:,:), m(:), t(:)
  integer               :: runtime
  
  allocate(S(L+2,L+2),m(n_meas),t(n_meas),BE(n_meas))
  
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S)
  call animate_lattice(S,'')
  
  call run_sim(S,BE,BJ,h,t,m,runtime)
  
  call close_lattice_plot()
  call results_out(BJ,BE(n_meas),h,runtime)
  call line_plot(real(t,dp),BE,'t','energy','','',1)
  call line_plot(real(t,dp),real(m,dp),'t','magnetization','','',2)
  
  deallocate(S,m,t,BE)
end program
