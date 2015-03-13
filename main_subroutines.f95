module main_subroutines
  use constants 
  implicit none
  private
  public :: calc_energy
contains
  subroutine calc_energy(BE,S,BJ,h)
    real(dp), intent(out) :: BE
    integer, intent(in) :: S(:,:)
    real(dp), intent(in) :: h, BJ
    integer :: i, j

    if (size(S,1) < 2) return !check
    
    BE = 0._dp ! initialze energy 
    
    do i = 2,L+1
      do j = 2,L+1
        BE = BE - BJ*S(i,j)*(S(i-1,j) + S(i+1,j) + S(i,j-1) + S(i,j+1))
      enddo
    enddo

    BE = 0.5_dp*BE ! account for double counting of pairs
    BE = BE - h*sum(S) ! add external field
  end subroutine
end module
