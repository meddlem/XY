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

  integer, parameter :: sweeps = 10000 ! number of sweeps through lattice 
  integer, parameter :: steps = sweeps*N ! number of iterations
end module
