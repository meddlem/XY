module proc_output
  use constants
  implicit none
  private
  public :: lin_fit, calc_chi, calc_h_mod

contains
  pure subroutine calc_chi(N,N_SWC,Chi,Chi_err)
    integer, intent(in)   :: N_SWC(:), N
    real(dp), intent(out) :: Chi, Chi_err

    Chi = sum(real(N_SWC,dp)/real(N,dp))/n_meas ! magnetic susc
    Chi_err = std_err(real(N_SWC,dp)/real(N,dp))
  end subroutine

  pure subroutine calc_h_mod(G,h_mod,h_mod_err)
    real(dp), intent(in)  :: G(:)
    real(dp), intent(out) :: h_mod, h_mod_err

    h_mod = sum(G)/n_meas ! helicity modulus 
    h_mod_err = std_err(G)
  end subroutine

  pure function std_err(A)
    ! calculates std error of A from blocked data
    real(dp), intent(in) :: A(:)
    real(dp) :: std_err, sigma_blocks_2, Avg(n_blocks)

    Avg = block_avg(A)
  
    sigma_blocks_2 = sum((Avg - sum(Avg)/n_blocks)**2)/(n_blocks-1)
    std_err = sqrt(sigma_blocks_2/n_blocks)
  end function

  pure function block_avg(A)
    ! returns array containing block average of A
    real(dp), intent(in) :: A(:)
    real(dp) :: block_avg(n_blocks)
    integer :: j

    do j = 0,(n_blocks-1)
      block_avg(j+1) = sum(A(n_avg*j+1:n_avg*(j+1)))/n_avg
    enddo
  end function

  pure subroutine lin_fit(slope,err_slope,offset,y,x)
    real(dp), intent(out) :: slope, err_slope, offset
    real(dp), intent(in) :: y(:), x(:)
    real(dp) :: mu_x, mu_y, ss_yy, ss_xx, ss_yx, s
    ! linear regression
    ! see also: http://mathworld.wolfram.com/LeastSquaresFitting.html

    mu_x = sum(x)/size(x)
    mu_y = sum(y)/size(y)

    ss_yy = sum((y - mu_y)**2)
    ss_xx = sum((x - mu_x)**2)
    ss_yx = sum((x - mu_x)*(y - mu_y))
    
    slope = ss_yx/ss_xx
    offset = mu_y - slope*mu_x
    s = sqrt((ss_yy - slope*ss_yx)/(size(x)-2))

    err_slope = s/(sqrt(ss_xx)*6._dp) ! error in slope calc
  end subroutine
end module
