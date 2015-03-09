program main
  use constants
  use initialize
  use markov
  use plotroutines
  use io
  implicit none
  
  ! variables:
  ! BE = beta*energy,
  ! dE = change in energy between new and initial config
  ! h = external field
  ! array containing Spins indexed as row,column

  real(dp) :: BJ, BE, BE_init, dE, h, dE_vals(9,2), BF_vals(9,2)
  integer, allocatable :: S(:,:)
  integer :: i, j, start_time, end_time, runtime
  
  allocate(S(L+2,L+2))
  call user_in(BJ,h)

  ! calculate values of dE and boltzmann factor 
  forall(i=1:9,j=1:2) dE_vals(i,j) = - 2._dp*BJ*(i-5) - 2._dp*h*(j*2-3) 
  BF_vals = exp(-dE_vals)
  
  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE,S,BJ,h)
  
  BE_init = BE
  call gnu_lattice_plot(S,1,'initial state')
 
  call system_clock(start_time)
  do i=1,100*N
    call gen_config(S,dE,dE_vals,BF_vals)
    BE = BE + dE
  enddo
  call system_clock(end_time)
  
  runtime = (end_time - start_time)/1000
  call results_out(BJ,BE,BE_init,h,runtime)
  call gnu_lattice_plot(S,2,'final state')
end program
