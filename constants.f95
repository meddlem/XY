module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)

  integer, parameter :: L = 100 ! lattice side
  integer, parameter :: N = L**2 ! number of spins

  integer, parameter :: meas_step = 100 ! interval between measurements
  integer, parameter :: steps = meas_step*N ! number of iterations
  integer, parameter :: n_meas = steps/meas_step ! total number of measurements
  integer, parameter :: meas_start = 100 ! start measurement after .. 

  integer, parameter :: n_corr = 30 ! number of spins to calc correlation over (diagonal elements)
  integer, parameter :: r_max = 10 ! distances over which to calc correlation function
end module
