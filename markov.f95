module markov 
  use constants
  implicit none
  private
  public :: gen_config

contains
  
  subroutine gen_config(S,dE,dE_vals,BF_vals)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(out) :: dE
    real(dp), intent(in) :: dE_vals(:,:), BF_vals(:,:)
    integer, allocatable :: S_t(:,:)
    integer :: i, j, ind_BJ, ind_h, x(2), S_nbrs
    real(dp) :: r, BF, dE_t
    
    allocate(S_t(L+2,L+2))
    S_t = S ! initialize trial config
    dE = 0._dp ! init dE

    ! create trial config by flipping 1 spin:
    call random_spin(x) 
    i = x(1)+1; j = x(2)+1 ! adjust for zero padding
    S_t(i,j) = -S_t(i,j) ! flip
    S_nbrs = S_t(i-1,j) + S_t(i+1,j) + S_t(i,j-1) + S_t(i,j+1) 
    
    ! get dE from matrix
    ind_BJ = S_t(i,j)*S_nbrs + 5
    ind_h = (S_t(i,j) + 3)/2
    dE_t = dE_vals(ind_BJ,ind_h)
   
    if (dE_t < 0._dp) then  
      S = S_t ! if energy decreases always accept
      dE = dE_t
    else ! else accept config with probability of BF
      BF = BF_vals(ind_BJ,ind_h)
      call random_number(r)
      if (r<BF) then
        S = S_t
        dE = dE_t
      endif
    endif
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
