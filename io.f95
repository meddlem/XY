module io
  use constants
  implicit none
  private
  public :: user_in, results_out
contains

  subroutine user_in(BJ)
    real(dp), intent(out) :: BJ
  
    write(*,'(/,A,/)') '************ Input *************' 
    write(*,'(A)',advance='no') "Beta*J = " 
    read(*,*) BJ
    write(*,'(A)') "Running simulation..."
  end subroutine

  subroutine results_out(BJ,BE,runtime) 
    real(dp), intent(in) :: BJ, BE
    integer, intent(in) :: runtime

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "Beta*J :", BJ
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Final Energy", BE
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    call system('cat output.txt')
  end subroutine
end module
