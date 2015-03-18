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

  real(dp), allocatable :: BE(:), c_ss(:), r(:)
  real(dp)              :: BJ, h
  integer, allocatable  :: S(:,:), m(:), t(:)
  integer               :: runtime
  
  allocate(S(L+2,L+2),m(n_meas),t(n_meas),BE(n_meas),c_ss(r_max),r(r_max))
  
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S)
  call animate_lattice('')
  
  call run_sim(S,BE,BJ,h,t,r,m,runtime,c_ss)
  
  call close_lattice_plot()
  call results_out(BJ,BE(n_meas),h,runtime)
  call line_plot(real(t,dp),BE,'t','energy','','',1)
  call line_plot(real(t,dp),real(m,dp),'t','magnetization','','',2)
  call line_plot(log(r),-log(c_ss),'log r','-log corr','','',3)
  
  deallocate(S,m,t,r,BE)
end program
