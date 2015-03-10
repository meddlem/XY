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
  ! dE_vals: contains all possible values of energy change
  ! BF_vals: all possible values of the boltzmann factor of dE

  real(dp) :: BJ, BE, BE_init, dE, h, dE_vals(9,2), BF_vals(9,2)
  integer, allocatable :: S(:,:)
  integer :: i, start_time, end_time, runtime
  
  allocate(S(L+2,L+2))
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE,S,BJ,h)
  call init_vals(dE_vals,BF_vals,BJ,h)

  BE_init = BE
  call lattice_plot(S,1,'',.true.)
 
  call system_clock(start_time)
  do i=1,steps
    call gen_config(S,dE,dE_vals,BF_vals)
    BE = BE + dE
    
    if (mod(i,5*N)==0) call write_lattice(S) !plot every sweep
  enddo
  call system_clock(end_time)

  call system('pkill gnuplot') !needed for now
  
  runtime = (end_time - start_time)/1000
  call lattice_plot(S,2,'final state',.false.)
  call results_out(BJ,BE,BE_init,h,runtime)
  deallocate(S)
end program
