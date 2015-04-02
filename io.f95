module io
  use constants
  use plotroutines
  implicit none
  private
  public :: get_usr_args, user_in, results_out
contains

  subroutine get_usr_args(use_T)
    logical, intent(out) :: use_T
    
    character(10) :: arg
    integer       :: i

    use_T = .false. ! by default use betaJ for input
    
    do i=1,iargc()
      call getarg(i,arg)
      if (trim(arg) == '-T') then
        use_T = .true.
      endif
    enddo
  end subroutine

  subroutine user_in(BJ,L,use_T)
    real(dp), intent(out) :: BJ
    integer, intent(out)  :: L
    logical, intent(in)   :: use_T

    real(dp) :: T
  
    write(*,'(/,A,/)') '************ Input *************' 
    if (use_T) then
      write(*,'(A)',advance='no') "T = " 
      read(*,*) T
      BJ = 1._dp/T
    else
      write(*,'(A)',advance='no') "BJ = " 
      read(*,*) BJ
    endif

    write(*,'(A)',advance='no') "L = " 
    read(*,*) L
    write(*,'(A)') "Running simulation..."

  end subroutine

  subroutine results_out(BJ,t,BE,h_mod,Chi,runtime) 
    real(dp), intent(in) :: BJ, t(:), BE(:), h_mod, Chi
    integer, intent(in) :: runtime

    character(30) :: rowfmt
    logical       :: exs
    
    write(rowfmt, '(A)') '(F7.5,3X,F7.5)'

    open(12,access = 'sequential',file = 'output.txt')
      write(12,'(/,A,/)') '*********** Summary ***********' 
      write(12,*) "BJ :", BJ
    
      write(12,'(/,A,/)') '*********** Output ************' 
      write(12,'(A,I6,A)') "Runtime : ", runtime, " s"
      write(12,*) "Helicity modulus", h_mod
      write(12,*) "Magnetic Susceptibility", Chi
      write(12,'(/,A,/)') '*******************************' 
    close(12)
    
    ! append spin stiffness calculation result to file
    inquire(file='TvsHmod.dat',exist=exs)
    if (exs) then
      open(12,file ='TvsHmod.dat',status='old',position='append',&
        action='write')
    else 
      open(12,file ='TvsHmod.dat',status='new',action='write')
    endif
      write(12,rowfmt) 1._dp/BJ, h_mod
    close(12)

    call line_plot(real(t,dp),BE,'t','energy','','',1)
    call system('cat output.txt')
  end subroutine
end module
