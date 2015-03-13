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

  real(dp) :: BJ, p, h, BE_tmp, m_tmp, dE
  real(dp), allocatable :: t(:), BE(:), m(:)
  integer, allocatable :: S(:,:)
  integer :: i, j, start_time, end_time, runtime
  
  allocate(S(L+2,L+2),m(sweeps+1),t(sweeps+1),BE(sweeps+1))
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S)
  call init_energy(BE_tmp,S,BJ,h)
  call animate_lattice(S,'')

  ! initialize some needed variables
  j = 1
  m(j) = sum(S)
  BE(j) = BE_tmp
  t = (/(i,i=0,sweeps)/)
  p = 1 - exp(-2._dp*BJ)

  call system_clock(start_time)
  do i=1,steps
    call gen_config(S,m_tmp,dE,p)
    BE_tmp = BE_tmp + dE

    if (mod(i,N)==0) then
      j = j+1
      m(j) = m_tmp
      BE(j) = BE_tmp ! record energy every sweep
    endif
    if (mod(i,N/10)==0) call write_lattice(S) ! write lattice to pipe
  enddo
  call system_clock(end_time)

  runtime = (end_time - start_time)/1000

  call close_lattice_plot()
  call results_out(BJ,BE_tmp,BE(1),h,runtime)
  call line_plot(t,m,'t','magnetization','','',1)
  deallocate(S)
end program
