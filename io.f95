module io
  use constants
  use plotroutines
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(J,L)
    real(dp), intent(out) :: J
    integer, intent(out)  :: L
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "J = " 
    read(*,*) J
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L
    write(*,'(A)') "Running simulation..."
  end subroutine

  subroutine results_out(J,t,E,h_mod,Chi,runtime) 
    real(dp), intent(in) :: J, t(:), E(:), h_mod, Chi
    integer, intent(in) :: runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "J :", J
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Helicity modulus", h_mod
      write(12,*) "Magnetic Susceptibility", Chi
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    call line_plot(real(t,dp),E,'t','energy','','',1)
    call system('cat output.txt')
  end subroutine
end module
