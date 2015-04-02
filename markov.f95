module markov 
  use constants
  use proc_output
  use plotroutines
  implicit none
  private
  public :: run_sim

contains
  subroutine run_sim(S,BE,BJ,t,h_mod,h_mod_err,Chi,Chi_err,runtime)
    real(dp), intent(inout) :: S(:,:,:), t(:), BE(:), BJ
    integer, intent(out)    :: runtime
    real(dp), intent(out)   :: Chi, Chi_err, h_mod, h_mod_err

    real(dp), allocatable :: G(:)
    integer, allocatable  :: N_SWC(:)
    integer   :: i, j, L, N, start_time, end_time, N_SWC_tmp
    
    allocate(G(n_meas),N_SWC(n_meas))
    ! initialize needed variables
    j = 0
    L = size(S,2)
    N = L**2
    t = real((/(i,i=0,n_meas-1)/),dp)
    
    call animate_lattice(S)
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

    runtime = (end_time - start_time)/1000 ! calculate runtime
    call calc_chi(N,N_SWC,Chi,Chi_err)
    call calc_h_mod(G,h_mod,h_mod_err)
    call line_plot(t,G,'','','','',2)
    deallocate(G,N_SWC)
  end subroutine

  subroutine gen_config(S,L,BJ,N_SWC)
    real(dp), intent(inout) :: S(:,:,:)
    integer, intent(out)    :: N_SWC
    integer, intent(in)     :: L
    real(dp), intent(in)    :: BJ

    integer, allocatable :: C(:,:)
    logical, allocatable :: in_cluster(:,:)
    integer   :: i, j, ss(2), nn(4,2)
    real(dp)  :: u(2)
    
    ! initialize variables 
    allocate(C(L**2,2),in_cluster(L**2,L**2))
    C = 0 ! cluster
    in_cluster = .false. ! tags spins 
    i = 1 ! labels for spin in cluster
    N_SWC = 1 ! number of spins in cluster

    call random_idx(ss,L) ! start cluster by choosing 1 spin
    call random_dir(u) ! get random unit vector

    C(1,:) = ss ! add chosen spin to cluster  
    in_cluster(ss(1),ss(2)) = .true. ! tag spin     
    call flip(S,ss(1),ss(2),u) ! flip initial spin
    
    do while (i<=N_SWC)
      ss = C(i,:) ! pick a spin x in the cluster
      nn = get_nn_idx(ss,L) ! get nearest neighbors of spin ss
      
      do j = 1,4 ! iterate over nearest neighbors of ss
        call try_add(S,C,in_cluster,N_SWC,u,ss,nn(j,:),BJ)
      enddo
      i = i+1 ! move to next spin in cluster
    enddo
    deallocate(C,in_cluster)
  end subroutine

  subroutine try_add(S,C,in_cluster,N_SWC,u,ss,nn,BJ)
    real(dp), intent(inout) :: S(:,:,:)
    integer, intent(inout)  :: N_SWC, C(:,:)
    logical, intent(inout)  :: in_cluster(:,:)
    integer, intent(in)     :: nn(:), ss(:)
    real(dp), intent(in)    :: BJ, u(:)

    integer   :: i, j, i_0, j_0
    real(dp)  :: r, p, S1_dot_u, S2_dot_u
    
    ! initial spin lattice coords
    i_0 = ss(1) 
    j_0 = ss(2)

    ! neigbor spin
    i = nn(1)
    j = nn(2)
    
    if (.not. in_cluster(i,j)) then ! check if spin already visited 
      call random_number(r)
      S1_dot_u = dot_product(S(:,i_0,j_0),u)
      S2_dot_u = dot_product(S(:,i,j),u)
      p = 1 - exp(2._dp*BJ*S1_dot_u*S2_dot_u)

      if (r < p) then ! add spin to cluster with probability p
        N_SWC = N_SWC+1 ! increase nr of spins in cluster
        
        ! add to cluster
        C(N_SWC,:) = [i,j] 
        in_cluster(i,j) = .true. 

        call flip(S,i,j,u) 
      endif
    endif
  end subroutine

  pure function get_nn_idx(x,L)
    ! returns indices of nearest neighbors of x_ij, accounting for PBC
    integer, intent(in) :: x(2), L
    integer :: get_nn_idx(4,2)

    get_nn_idx(1,:) = merge(x + [1,0], [1,x(2)], x(1) /= L)
    get_nn_idx(2,:) = merge(x + [0,1], [x(1),1], x(2) /= L) 
    get_nn_idx(3,:) = merge(x - [1,0], [L,x(2)], x(1) /= 1) 
    get_nn_idx(4,:) = merge(x - [0,1], [x(1),L], x(2) /= 1) 
  end function

  pure subroutine flip(S,i,j,u)
    ! flips spin S(i,j) wrt unit vector u
    real(dp), intent(inout) :: S(:,:,:)
    real(dp), intent(in)    :: u(:)
    integer, intent(in)     :: i, j

    real(dp) :: S_tmp(2)
    
    S_tmp = S(:,i,j) 
    
    S_tmp = S_tmp - 2._dp*dot_product(S_tmp,u)*u  
    S_tmp = S_tmp/sqrt(sum(S_tmp**2)) ! ensure normalization of spins   
    
    S(:,i,j) = S_tmp
  end subroutine
  
  subroutine random_idx(x,L)
    ! returns index of randomly picked spin
    integer, intent(in)  :: L
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = (L-1)*u + 1
    x = int(u) 
  end subroutine

  subroutine random_dir(r)
    ! returns random unit vector 
    real(dp), intent(out) :: r(:)
    real(dp) :: u

    call random_number(u)
    u = 2._dp*pi*u
    r = [cos(u),sin(u)]
  end subroutine
  
  pure subroutine calc_energy(BE,S,L,BJ)
    ! calculate energy of the system
    real(dp), intent(out) :: BE
    real(dp), intent(in)  :: S(:,:,:), BJ
    integer, intent(in)   :: L
    integer               :: i, j, k, nn(4,2)
    
    BE = 0._dp ! init energy 

    do i = 1,L
      do j = 1,L
        nn = get_nn_idx([i,j],L) 
        do k = 1,4
          BE = BE - BJ*dot_product(S(:,i,j),S(:,nn(k,1),nn(k,2)))
        enddo
      enddo
    enddo

    BE = 0.5_dp*BE ! correct for double counting of pairs
  end subroutine

  pure subroutine helicity_mod(G,S,L,BE,BJ)
    ! calculates helicity modulus/spin wave stifness
    real(dp), intent(out) :: G
    real(dp), intent(in)  :: S(:,:,:), BJ, BE
    integer, intent(in)   :: L
    
    real(dp), allocatable :: dtheta_x(:,:), dtheta_y(:,:)
    integer               :: i, j, N
    
    ! initialize variables
    allocate(dtheta_x(L,L),dtheta_y(L,L))
    G = 0._dp
    N = L**2
    
    ! calculate angle difference between neighboring spins in x and y dirs
    do i=1,L
      do j=1,L
        dtheta_x(i,j) = angle(S(:,i,j)) - angle(S(:,modulo(i,L)+1,j))
        dtheta_y(i,j) = angle(S(:,i,j)) - angle(S(:,i,modulo(j,L)+1))
      enddo
    enddo

    G = (-BE/BJ - BJ*(sum(sin(dtheta_x))**2 + sum(sin(dtheta_y))**2))/(2._dp*N)
    deallocate(dtheta_x,dtheta_y)
  end subroutine
  
  pure function angle(S)
    real(dp), intent(in) :: S(:)
    real(dp)             :: angle
    ! convert spin vector to angle wrt x-axis
    angle = atan2(S(2),S(1)) 
  end function
end module
