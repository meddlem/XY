module io
  use constants
  use plotroutines
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(BJ,L)
    real(dp), intent(out) :: BJ
    integer, intent(out)  :: L
    real(dp) :: T
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "T = " 
    read(*,*) T
    write(*,'(A)',advance='no') "L = " 
    read(*,*) L
    write(*,'(A)') "Running simulation..."

    BJ = 1._dp/T
  end subroutine

  subroutine results_out(BJ,t,BE,h_mod,Chi,runtime) 
    real(dp), intent(in) :: BJ, t(:), BE(:), h_mod, Chi
    integer, intent(in) :: runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "BJ :", BJ
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Helicity modulus", h_mod
      write(12,*) "Magnetic Susceptibility", Chi
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    call line_plot(real(t,dp),BE,'t','energy','','',1)
    call system('cat output.txt')
  end subroutine
end module
