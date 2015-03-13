module markov 
  use constants
  implicit none
  private
  public :: gen_config

contains
  subroutine gen_config(S,m,dE,p)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(out) :: m, dE
    real(dp), intent(in) :: p
    integer, allocatable :: C(:,:)
    integer :: i, j, S_init, s_cl, x(2), nn(4,2)
    real(dp) :: r
    
    allocate(C(N,2))

    ! initialize variables 
    dE = 0._dp ! init dE, 
    !need to implement calculation of dE based on bonds with neighbors?
    i = 1 ! labels spin in cluster
    s_cl = 1 ! number of spins in cluster
    C = 0 ! init array that holds indices of all spins in cluster

    call random_spin(x) ! start cluster by choosing 1 spin

    S_init = S(x(1),x(2)) ! save state of chosen spin
    C(1,:) = x ! add spin to cluster     
    S(x(1),x(2)) = -S_init ! flip initial spin
    
    do while (i<=s_cl)
      x = C(i,:) ! pick a spin x in the cluster
      nn = nn_idx(x) ! get nearest neighbors of spin x
      
      ! iterate over neighbors of x
      do j = 1,4 
        if (S(nn(j,1),nn(j,2))==S_init) then 
          call random_number(r)

          if (r<p) then ! add spin to cluster with probability p
            s_cl = s_cl+1
            C(s_cl,:) = nn(j,:) 
            
            S(nn(j,1),nn(j,2)) = -S_init ! flip spin
          endif
        endif
      enddo
      i = i+1 ! move to next spin in cluster
    enddo 

    m = sum(S) ! calculate instantaneous magnetization
    deallocate(C)
  end subroutine

  subroutine random_spin(x)
    ! returns index of randomly picked spin
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = L*u + 0.5_dp
    x = nint(u) ! index of spin to flip
    x = x + 1 ! adjust for zero padding
  end subroutine

  function nn_idx(x)
    ! returns indices of nearest neighbors of x_ij
    integer, intent(in) :: x(2)
    integer :: nn_idx(4,2)

    nn_idx(1,:) = x + [1,0] 
    nn_idx(2,:) = x + [0,1] 
    nn_idx(3,:) = x - [1,0] 
    nn_idx(4,:) = x - [0,1] 
  end function
end module
