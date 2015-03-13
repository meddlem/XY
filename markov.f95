module markov 
  use constants
  use plotroutines
  implicit none
  private
  public :: run_sim, gen_config

contains
  subroutine run_sim(S,BE,BJ,h,t,m,runtime)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(inout) :: BE(:) 
    integer, intent(out) :: t(:), m(:), runtime
    real(dp), intent(inout) :: BJ, h
    real(dp) :: p
    integer :: i, j, start_time, m_tmp, end_time
  
    ! initialize needed variables
    j = 0
    h = 0._dp ! overwrite user setting, just in case 
    t = (/(i,i=0,n_meas-1)/)
    p = 1 - exp(-2._dp*BJ)

    call system_clock(start_time)
    do i=1,steps
      call gen_config(S,m_tmp,p)

      if (mod(i,meas_step)==0) then
        j = j+1
        m(j) = m_tmp
        call calc_energy(BE(j),S,BJ,h)
      endif
      if (mod(i,N/10)==0) call write_lattice(S) ! write lattice to pipe
    enddo
    
    call system_clock(end_time)
    runtime = (end_time - start_time)/1000
  end subroutine

  subroutine gen_config(S,m,p)
    integer, intent(inout) :: S(:,:)
    integer, intent(out) :: m
    real(dp), intent(in) :: p
    integer, allocatable :: C(:,:)
    integer :: i, j, S_init, s_cl, x(2), nn(4,2)
    real(dp) :: r
    
    allocate(C(N,2))

    ! initialize variables 
    i = 1 ! labels spin in cluster
    s_cl = 1 ! number of spins in cluster
    C = 0 ! init array that holds indices of all spins in cluster

    call random_spin(x) ! start cluster by choosing 1 spin

    S_init = S(x(1),x(2)) ! save state of chosen spin
    C(1,:) = x ! add chosen spin to cluster     
    S(x(1),x(2)) = -S_init ! flip initial spin
    
    do while (i<=s_cl)
      x = C(i,:) ! pick a spin x in the cluster
      nn = nn_idx(x) ! get nearest neighbors of spin x
      
      ! iterate over neighbors of x
      do j = 1,4 
        if (S(nn(j,1),nn(j,2)) == S_init) then 
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
