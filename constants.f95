module constants
  implicit none
  ! module contains all constants used in the program
  
  ! dp: compiler specific kind for double precision float
  ! lng: compiler specific kind for long integer

  ! NOTE: IF YOU MAKE ANY CHANGES HERE RECOMPILE ALL MODULES: "make -B" 
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: lng = selected_int_kind(8)
  real(dp), parameter :: pi = 4._dp*atan(1._dp)

  integer, parameter :: steps = 20000 
  integer, parameter :: meas_start = 5000 
  integer, parameter :: n_meas = steps - meas_start 
  integer, parameter :: plot_interval = 99 
end module
