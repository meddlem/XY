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
    integer, allocatable :: S0(:,:), S_t(:,:)
    integer :: i, j, x(2), S_nbrs
    real(dp) :: BF, dE_t, r
    
    allocate(S0(L+2,L+2), S_t(L,L))
    
    S_t = S ! initialize trial config
    dE = 0._dp ! init dE
    call random_spin(x)
    
    ! create trial config by flipping 1 spin:
    i = x(1); j = x(2)
    S_t(i,j) = -S_t(i,j)
    
    ! add zero padding, probably best to do outside of this module?
    S0 = 0
    S0(2:L+1,2:L+1) = S_t
    i = i+1; j = j+1 ! adjust indices accordingly 
    
    ! calculate change in BE
    S_nbrs = S0(i-1,j) + S0(i+1,j) + S0(i,j-1) + S0(i,j+1)
    dE_t = -2._dp*BJ*S0(i,j)*S_nbrs - 2._dp*h*S0(i,j) 
   
    if (dE_t < 0._dp) then  
      S = S_t ! if energy decreases always accept
      dE = dE_t
    else ! else accept config with probability of BF
      BF = exp(-dE_t)
      call random_number(r)
      if (r<BF) then
        S = S_t
        dE = dE_t
      endif
    endif

    deallocate(S0,S_t)
  end subroutine

  subroutine random_spin(x)
    ! returns index of randomly picked spin
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = L*u + 0.5_dp
    x = nint(u) ! index of spin to flip
  end subroutine
end module
