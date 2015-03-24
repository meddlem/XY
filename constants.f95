module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)
  real(dp), parameter :: pi = 4._dp*atan(1._dp)

  integer, parameter :: L = 32 ! lattice side
  integer, parameter :: N = L**2 ! number of spins

  integer, parameter :: meas_step = 100 ! interval between measurements
  integer, parameter :: steps = 100000 ! number of iterations
  integer, parameter :: n_meas = steps/meas_step ! total number of measurements
  integer, parameter :: meas_start = 1000 ! start measurement after .. steps 
  integer, parameter :: plot_interval = 100 ! plot every .. steps
end module
