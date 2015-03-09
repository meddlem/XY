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
  
  allocate(S(L+2,L+2))
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE,S,BJ,h)
  
  BE_init = BE
  call lattice_plot(S,1,'',.true.)
 
  call system_clock(start_time)
  do i=1,steps
    call gen_config(S,dE,BJ,h)
    if (mod(i,N)==0) then
      call write_lattice(S) 
    endif
    BE = BE + dE
  enddo
  call system_clock(end_time)
  
  call system('pkill gnuplot')
  
  runtime = (end_time - start_time)/1000
  call results_out(BJ,BE,BE_init,h,runtime)
  call lattice_plot(S,1,'final state',.false.)
end program
