module markov 
  use constants
  implicit none
  private
  public :: gen_config

contains
  
  subroutine gen_config(S,dE,BJ,h)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(out) :: dE
    real(dp), intent(in) :: BJ, h
    integer :: i, j, S0(L+2,L+2), S_t(L,L)
    real(dp) :: x(2), BF, dE_t, r
    
    S_t = S

    ! create trial config by flipping 1 spin:
    call random_number(x)
    x = L*x + 0.5_dp
    i = nint(x(1)) ! index of spin to flip
    j = nint(x(2)) ! index of spin to flip
    print *, i, j

    S_t(i,j) = -S_t(i,j) ! flip spin
    
    S0 = 0
    S0(2:L+1,2:L+1) = S_t ! add zero padding
    i = i+1
    j = j+1 !adjust i, j accordingly

    ! you can just save the possible energy changes, boltzmann factors in &
    ! an array, may save some computation time. could do this in init_energy
    dE_t = -2._dp*BJ*S0(i,j)*(S0(i-1,j) + S0(i+1,j) + S0(i,j-1) + S0(i,j+1)) - 2._dp*h*S0(i,j) ! calculate change in energy
    
    ! accept trial config 
    if (dE_t < 0._dp) then 
      S = S_t ! if energy decreases always accept
      dE = dE_t
    else ! else accept config with probability of BF
      BF = exp(-dE)
      call random_number(r)
      
      if (r<BF) then
        S = S_t
        dE = dE_t
      else
        S = S
        dE = 0._dp
      endif
    endif
  end subroutine
end module
