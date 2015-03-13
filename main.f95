program main
  use constants
  use initialize
  use main_subroutines
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

  real(dp) :: BJ, p, h, m_tmp
  real(dp), allocatable :: t(:), BE(:), m(:)
  integer, allocatable :: S(:,:)
  integer :: i, j, start_time, end_time, runtime
  
  allocate(S(L+2,L+2),m(n_meas+1),t(n_meas+1),BE(n_meas+1))
  call user_in(BJ,h)
  call init_random_seed()
  call init_lattice(S)
  call calc_energy(BE(1),S,BJ,h)
  call animate_lattice(S,'')

  ! initialize some needed variables
  j = 1
  h = 0._dp ! overwrite user setting, just in case 
  m(j) = sum(S)
  t = (/(i,i=0,n_meas)/)
  p = 1 - exp(-2._dp*BJ)

  call system_clock(start_time)
  do i=1,steps
    call gen_config(S,m_tmp,p)

    if (mod(i,meas_step)==0) then
      j = j+1
      m(j) = m_tmp
      call calc_energy(BE(j),S,BJ,h)
    endif
    if (mod(i,N/10)==0) call write_lattice(S) ! write lattice to pipe
  enddo
  call system_clock(end_time)

  runtime = (end_time - start_time)/1000

  call close_lattice_plot()
  !call results_out(BJ,BE_tmp,BE(1),h,runtime)
  call line_plot(t,BE,'t','energy','','',1)
  call line_plot(t,m,'t','magnetization','','',2)
  deallocate(S)
end program
