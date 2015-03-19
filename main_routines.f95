module main_routines
  use constants
  implicit none
  private 
  public :: s_corr, lin_fit

contains

  pure subroutine s_corr(g,S)
    real(dp), intent(out) :: g(:)
    integer, intent(in) :: S(:,:)
    real(dp) :: g_tmp(n_corr,r_max)
    integer :: i, r_0, r_1
     
    r_0 = (L-n_corr)/2
    r_1 = r_0 + 1

    do i=1,n_corr
      g_tmp(i,:) = S(i+r_0,i+r_0)*&
        (S(i+r_0,i+r_1:i+r_0+r_max) + S(i+r_1:i+r_0+r_max,i+r_0))/2._dp
    enddo

    g = sum(g_tmp,1)/n_corr 
  end subroutine

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
