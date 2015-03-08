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

  real(dp) :: BJ, BE, BE_init, dE, h
  integer, allocatable :: S(:,:)
  integer :: i, start_time, end_time, runtime
  
  allocate(S(L,L))
  call user_in(BJ,h)
  ! dE_vals = 2._dp*BJ*(/(i,i=-4,4)/) ! possible values of dE
  ! BF_vals = exp(-dE_vals) ! possible values of boltzmann factor

  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE,S,BJ,h)
  
  BE_init = BE
  call gnu_lattice_plot(S,1,'initial state')
 
  call system_clock(start_time)
  do i=1,100*N
    call gen_config(S,dE,BJ,h)
    BE = BE + dE
  enddo
  call system_clock(end_time)
  
  runtime = (end_time - start_time)/1000
  call results_out(BJ,BE,BE_init,h,runtime)
  call gnu_lattice_plot(S,2,'final state')
end program
