module markov 
  use constants
  use main_routines
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(S,BE,BJ,t,h_mod,Chi,runtime)
    real(dp), intent(inout) :: S(:,:,:), t(:), BE(:), BJ
    integer, intent(out)    :: runtime
    real(dp), intent(out)   :: Chi

    real(dp), allocatable :: G(:)
    integer, allocatable :: N_SWC(:)
    real(dp)  :: h_mod
    integer   :: i, j, L, N, start_time, end_time, N_SWC_tmp
    
    allocate(G(n_meas),N_SWC(n_meas))
    ! initialize needed variables
    j = 0
    L = size(S,2)
    N = L**2
    t = real((/(i,i=0,n_meas-1)/),dp)
    
    call animate_lattice(L)
    call system_clock(start_time)
    do i=1,steps
      call gen_config(S,L,BJ,N_SWC_tmp)

      if (i>meas_start) then
        j = j+1
        N_SWC(j) = N_SWC_tmp
        call calc_energy(BE(j),S,L,BJ)
        call helicity_mod(G(j),S,L,BE(j),BJ)
      endif
  
      if (mod(i,plot_interval) == 0) then
        call write_lattice(S,L) ! lattice to pipe
      endif
    enddo    
    call system_clock(end_time)
    call close_lattice_plot()

    ! calculate runtime
    runtime = (end_time - start_time)/1000
    
    Chi = sum(real(N_SWC,dp)**2/real(N,dp)**2)/n_meas ! magnetic susc
    h_mod = sum(G)/n_meas ! helicity modulus 
    deallocate(G,N_SWC)
  end subroutine

  subroutine gen_config(S,L,BJ,N_SWC)
    real(dp), intent(inout) :: S(:,:,:)
    integer, intent(out)    :: N_SWC
    integer, intent(in)     :: L
    real(dp), intent(in)    :: BJ

    integer, allocatable :: C(:,:)
    logical, allocatable :: in_cluster(:,:)
    integer   :: i, j, x(2), nn(4,2)
    real(dp)  :: S_dot_u, u(2)
    
    ! initialize variables 
    allocate(C(L**2,2),in_cluster(L**2,L**2))
    C = 0 ! cluster
    in_cluster = .false. ! tags spins 
    i = 1 ! labels for spin in cluster
    N_SWC = 1 ! number of spins in cluster

    call random_idx(x,L) ! start cluster by choosing 1 spin
    call random_dir(u) ! get random unit vector

    C(1,:) = x ! add chosen spin to cluster  
    in_cluster(x(1),x(2)) = .true. ! tag spin     

    S(:,x(1),x(2)) = Flip(S(:,x(1),x(2)),u) ! flip initial spin
    S_dot_u = dot_product(S(:,x(1),x(2)),u) 
    
    do while (i<=N_SWC)
      x = C(i,:) ! pick a spin x in the cluster
      nn = nn_idx(x,L) ! get nearest neighbors of spin x
      
      do j = 1,4 ! iterate over nearest neighbors of x
        call try_add(S,C,in_cluster,N_SWC,S_dot_u,u,nn(j,:),BJ)
      enddo
      i = i+1 ! move to next spin in cluster
    enddo
    deallocate(C,in_cluster)
  end subroutine

  subroutine try_add(S,C,in_cluster,N_SWC,S_dot_u,u,s_idx,BJ)
    real(dp), intent(inout) :: S(:,:,:)
    integer, intent(inout)  :: N_SWC, C(:,:)
    logical, intent(inout)  :: in_cluster(:,:)
    integer, intent(in)     :: s_idx(:)
    real(dp), intent(in)    :: S_dot_u, BJ, u(:)

    integer   :: i, j 
    real(dp)  :: r, p, Sy_dot_u

    i = s_idx(1)
    j = s_idx(2)
    
    if (.not. in_cluster(i,j)) then ! check if spin already visited 
      call random_number(r)
      Sy_dot_u = dot_product(S(:,i,j),u)
      p = 1 - exp(min(2._dp*BJ*S_dot_u*Sy_dot_u,0._dp))

      if (r<p) then ! add spin to cluster with probability p
        N_SWC = N_SWC+1 ! increase nr of spins in cluster
        C(N_SWC,:) = s_idx ! add to cluster
        in_cluster(i,j) = .true. ! only check each bond once or multiple times?

        S(:,i,j) = Flip(S(:,i,j),u) ! flip spin 
      endif
    endif
  end subroutine

  pure function nn_idx(x,L)
    ! returns indices of nearest neighbors of x_ij, accounting for PBC
    integer, intent(in) :: x(2), L
    integer :: nn_idx(4,2)

    nn_idx(1,:) = merge(x + [1,0], [1,x(2)], x(1) /= L)
    nn_idx(2,:) = merge(x + [0,1], [x(1),1], x(2) /= L) 
    nn_idx(3,:) = merge(x - [1,0], [L,x(2)], x(1) /= 1) 
    nn_idx(4,:) = merge(x - [0,1], [x(1),L], x(2) /= 1) 
  end function

  pure function Flip(S,u)
    ! flips spin wrt unit vector u
    real(dp), intent(in) :: S(:), u(:)
    real(dp) :: Flip(2)
    
    Flip = S - 2._dp*dot_product(S,u)*u  
    Flip = Flip/sqrt(sum(Flip**2)) ! ensure normalization of spins   
  end function
  
  subroutine random_idx(x,L)
    ! returns index of randomly picked spin
    integer, intent(in)  :: L
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = (L-1)*u + 1
    x = int(u) ! index of spin to flip
  end subroutine

  subroutine random_dir(r)
    ! returns random unit vector 
    real(dp), intent(out) :: r(:)
    real(dp) :: u

    call random_number(u)
    u = 2._dp*pi*u
    r = [cos(u),sin(u)]
  end subroutine
  
  pure function angle(S)
    real(dp), intent(in) :: S(:)
    real(dp)             :: angle
    ! convert spin vector to angle wrt x-axis
    
    angle = atan2(S(2),S(1)) 
  end function

  pure subroutine calc_energy(BE,S,L,BJ)
    ! calculate energy of the system
    real(dp), intent(out) :: BE
    real(dp), intent(in)  :: S(:,:,:), BJ
    integer, intent(in)   :: L
    integer               :: i, j, k, nn(4,2)
    
    BE = 0._dp ! init energy 

    do i = 1,L
      do j = 1,L
        nn = nn_idx([i,j],L) ! get nearest neighbors of spin i,j
        do k = 1,4
          BE = BE - BJ*dot_product(S(:,i,j),S(:,nn(k,1),nn(k,2)))
        enddo
      enddo
    enddo

    BE = 0.5_dp*BE ! correct for double counting of pairs
  end subroutine

  pure subroutine helicity_mod(G,S,L,BE,BJ)
    real(dp), intent(out) :: G
    real(dp), intent(in)  :: S(:,:,:), BJ, BE
    integer, intent(in)   :: L
    
    real(dp), allocatable :: dthetax(:,:), dthetay(:,:)
    integer               :: i, j

    allocate(dthetax(L,L),dthetay(L,L))
    G = 0._dp

    do i=1,L
      do j=1,L
        dthetax(i,j) = angle(S(:,i,j)) - angle(S(:,modulo(i,L)+1,j))
        dthetay(i,j) = angle(S(:,i,j)) - angle(S(:,i,modulo(j,L)+1))
      enddo
    enddo

    G = -BE/(2._dp*BJ*L**2) ! double counting correction
    G = G - BJ*(sum(sin(dthetax))**2 + sum(sin(dthetay))**2)/(2._dp*L**2)
    deallocate(dthetax,dthetay)
  end subroutine
end module
