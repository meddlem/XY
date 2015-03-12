module markov 
  use constants
  implicit none
  private
  public :: gen_config

contains
  subroutine gen_config(S,dE,p)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(out) :: dE
    real(dp), intent(in) :: p
    integer, allocatable :: C_idx(:,:)
    integer :: i, j, S_tmp, s_cl, s_add, x(2), x_nn(2), nn(4,2)
    real(dp) :: r
    
    allocate(C_idx(N,2))

    ! initialize variables 
    dE = 0._dp ! init dE
    i = 1 ! labels spin in cluster
    s_cl = 1 ! number of spins in cluster
    s_add = 1 ! number of spins added to cluster
    C_idx = 0 ! init array that holds indices of all spins in cluster

    call random_spin(x) ! start cluster by choosing 1 spin
    C_idx(1,:) = x ! add spin to cluster     
    S_tmp = S(x(1),x(2)) ! save state of orig spin
    S(x(1),x(2)) = -S_tmp ! flip initial spin
    
    do while ((s_add /= 0) .or. (i<=s_cl))
      s_add = 0
      x = C_idx(i,:) ! pick a spin x in the cluster
      nn = nn_idx(x) ! get nearest neighbors of spin x
      
      ! iterate over neighbors of x
      do j = 1,4 
        x_nn = nn(j,:)
        if (S(x_nn(1),x_nn(2))==S_tmp) then 
          call random_number(r)

          if (r<p) then  
            s_cl = s_cl+1
            s_add = s_add+1
            C_idx(s_cl,:) = x_nn ! add spin to cluster with probability p

            S(x_nn(1),x_nn(2)) = -S_tmp ! flip spin so it's not visited again
          endif
        endif
      enddo
      i = i+1 ! move to next spin in cluster
    enddo 
    deallocate(C_idx)
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
