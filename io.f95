module io
  use constants
  use plotroutines
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(BK,L)
    real(dp), intent(out) :: BK
    integer, intent(out)  :: L
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "BK = " 
    read(*,*) BK
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L
    write(*,'(A)') "Running simulation..."
  end subroutine

  subroutine results_out(BK,t,BE,h_mod,Xi,runtime) 
    real(dp), intent(in) :: BK, t(:), BE(:), h_mod, Xi
    integer, intent(in) :: runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "BK :", BK
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Helicity modulus", h_mod
      write(12,*) "Magnetic Susceptibility", Xi
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    call line_plot(real(t,dp),BE,'t','energy','','',1)
    call system('cat output.txt')
  end subroutine
end module
