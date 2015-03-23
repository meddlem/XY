module markov 
  use constants
  use main_routines
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(S,BE,BK,t,m,runtime)
    real(dp), intent(inout) :: S(:,:,:), BE(:), BK, m(:)
    integer, intent(inout) :: t(:)
    integer, intent(out) :: runtime

    integer :: i, j, start_time, end_time
    real(dp) :: m_tmp
    
    ! initialize needed variables
    j = 0
    t = (/(i,i=0,n_meas-1)/)

    call system_clock(start_time)
    do i=1,steps
      call gen_config(S,m_tmp,BK)

      if (mod(i,meas_step) == 0) then
        j = j+1
        m(j) = m_tmp
        call calc_energy(BE(j),S,BK)
      endif

      if (mod(i,plot_interval) == 0) call write_lattice(S) ! write lattice to pipe
    enddo    
    call system_clock(end_time)
    runtime = (end_time - start_time)/1000
  end subroutine

  subroutine gen_config(S,m,BK)
    real(dp), intent(inout) :: S(:,:,:)
    real(dp), intent(in) :: BK
    real(dp), intent(out) :: m

    integer, allocatable :: C(:,:)
    logical, allocatable :: C_added(:,:)
    integer :: i, j, s_cl, x(2), nn(4,2)
    real(dp) :: sigma_x, u(2)
    
    allocate(C(N,2),C_added(N,N))
    ! initialize variables 
    i = 1 ! labels spin in cluster
    s_cl = 1 ! number of spins in cluster
    C_added = .false. ! tells us if spin was already considered for cluster
    C = 0 ! init array that holds indices of all spins in cluster

    call random_idx(x) ! start cluster by choosing 1 spin
    call random_dir(u) 
    C(1,:) = x
    C_added(x(1),x(2)) = .true. ! add chosen spin to cluster     
    
    ! flip initial spin
    sigma_x = dot_product(S(:,x(1),x(2)),u) 
    S(:,x(1),x(2)) = S(:,x(1),x(2)) - 2._dp*sigma_x*u  
    !sigma_x = dot_product(S(:,x(1),x(2)),u) 
    
    do while (i<=s_cl)
      x = C(i,:) ! pick a spin x in the cluster
      nn = nn_idx(x) ! get nearest neighbors of spin x
      
      do j = 1,4 ! iterate over neighbors of x
        call try_add(S,C,C_added,s_cl,sigma_x,u,nn(j,:),BK)
      enddo
      i = i+1 ! move to next spin in cluster
    enddo

    m = sum(S) ! calculate instantaneous magnetization
    deallocate(C,C_added)
  end subroutine

  subroutine try_add(S,C,C_added,s_cl,sigma_x,u,s_idx,BK)
    real(dp), intent(inout) :: S(:,:,:)
    integer, intent(inout) :: s_cl, C(:,:)
    logical, intent(inout) :: C_added(:,:)
    integer, intent(in) :: s_idx(:)
    real(dp), intent(in) :: sigma_x, BK, u(:)

    integer :: i, j 
    real(dp) :: r, p, b

    i = s_idx(1)
    j = s_idx(2)
    
    if (C_added(i,j) .eqv. .false.) then
      
      b = sigma_x*dot_product(S(:,i,j),u)
      p = 1 - exp(2*min(BK*b,0._dp)) ! check of dit echt klopt 
      call random_number(r)

      if (r<p) then ! add spin to cluster with probability p
        s_cl = s_cl+1 ! increase nr of spins in cluster
        C(s_cl,:) = s_idx ! add to cluster
        C_added(i,j) = .true. ! tag spin 

        S(:,i,j) = S(:,i,j) - 2._dp*dot_product(S(:,i,j),u)*u ! flip spin 
      endif
    endif
  end subroutine

  pure subroutine calc_energy(BE,S,BK)
    real(dp), intent(out) :: BE
    real(dp), intent(in) :: S(:,:,:), BK

    integer :: i, j, k, nn(4,2)
    
    BE = 0._dp ! initialze energy 

    do i = 1,L
      do j = 1,L
        nn = nn_idx([i,j]) ! get nearest neighbors of spin i,j
        do k = 1,4
          BE = BE - BK*dot_product(S(:,i,j),S(:,nn(k,1),nn(k,2)))
        enddo
      enddo
    enddo

    BE = 0.5_dp*BE ! account for double counting of pairs
  end subroutine
  
  pure function nn_idx(x)
    ! returns indices of nearest neighbors of x_ij, accounting for PBC
    integer, intent(in) :: x(2)
    integer :: nn_idx(4,2)

    nn_idx(1,:) = merge(x + [1,0], [1,x(2)], x(1) /= L)
    nn_idx(2,:) = merge(x + [0,1], [x(1),1], x(2) /= L) 
    nn_idx(3,:) = merge(x - [1,0], [L,x(2)], x(1) /= 1) 
    nn_idx(4,:) = merge(x - [0,1], [x(1),L], x(2) /= 1) 
  end function
  
  subroutine random_idx(x)
    ! returns index of randomly picked spin
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = L*u + 0.5_dp
    x = nint(u) ! index of spin to flip
  end subroutine

  subroutine random_dir(r)
    ! returns random unit vector 
    real(dp), intent(out) :: r(:)
    real(dp) :: u

    call random_number(u)
    u = 2._dp*pi*u
    r = [cos(u), sin(u)]
  end subroutine
end module
