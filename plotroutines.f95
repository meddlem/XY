module plotroutines
  use constants
  implicit none
  private
  public :: gnu_line_plot, gnu_lattice_plot 

contains
  subroutine gnu_line_plot(x,y1,xlabel,ylabel,label1,title,plot_no,y2,label2)
    real(dp), intent(in) :: x(:), y1(:)
    real(dp), intent(in), optional :: y2(:)
    character(*), intent(in) :: xlabel, ylabel, label1, title
    character(*), intent(in), optional :: label2
    integer, intent(in) :: plot_no
    character(1024) :: filename
    integer :: i, ret, m
    real(dp) :: xmin, xmax, xrange(2), ymin, ymax, yrange(2)
    
    if (size(y1)/=size(x)) print *, "error, arguments must be same size"
    
    if (present(y2)) then
      if (size(y1)/=size(y2)) then 
        print *, "error, arguments must be same size"
        return
      endif
    endif
    
    m = size(x)
    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y1)
    ymax = maxval(y1)

    xrange = [0._dp, xmax+(xmax-xmin)*0.1_dp]
    yrange = [ymin-(ymax-ymin)*0.1_dp, ymax+(ymax-ymin)*0.1_dp]

    open(10,access = 'sequential',file = 'xydata.dat')
    
    do i=1,m
      if (present(y2)) then
        write(10,*) x(i),y1(i),y2(i) ! write datapoints to file
      else
        write(10,*) x(i),y1(i) ! write datapoints to file
      endif
    enddo
    
    close(10,status = 'keep')
    
    ! create gnuplot command file
    write(filename,'(A,I1,A)') 'set output "plot',plot_no,'.png"'
    open(10,access = 'sequential',file = 'gplot.txt')
    
    ! set output terminal  
    write(10,*) 'set term pngcairo size 640,480 enhanced font "Verdana,10"'
    ! write(10,*) 'set term epscairo size 13cm,9cm font "Verdana,15"'
    write(10,*) filename
    ! set line color definitions
    write(10,*) &
      'set style line 1 lt 1 lc rgb "#ff0000" lw 2 #red'
    write(10,*) &
      'set style line 2 lt 1 lc rgb "#0000ff" lw 2 #blue'
    ! axes 
    write(10,*) 'set style line 11 lc rgb "#808080" lt 1'
    write(10,*) 'set border 31 back ls 11'
    write(10,*) 'set tics nomirror scale 0.75'
    write(10,*) 'set key right center'
    write(10,*) 'set mxtics 2'
    write(10,*) 'set mytics 2'
    ! grid 
    write(10,*) 'set style line 12 lc rgb "#808080" lt 0 lw 1'
    write(10,*) 'set grid back ls 12'
    ! plotrange
    write(10,*) 'set xrange [',xrange(1),':',xrange(2),']'
    write(10,*) 'set yrange [',yrange(1),':',yrange(2),']'
    ! plot labels
    write(10,*) 'set title "'//TRIM(title)//'"'
    write(10,*) 'set xlabel '//'"'//TRIM(xlabel)//'"'
    write(10,*) 'set ylabel '//'"'//TRIM(ylabel)//'"'
    
    if (m>0) then
      write(10,*) 'plot "xydata.dat" using 1:2 with line ls 10 t "", \'
      write(10,*) &
      ' "xydata.dat" using 1:2 with line ls 1 t "'//TRIM(label1)//'", \'
      if (present(y2)) then
        write(10,*) &
        ' "xydata.dat" using 1:3 with line ls 2 t "'//TRIM(label2)//'"'
      endif
    endif
    
    close(10,status = 'keep')

    ! now call gnuplot and plot the curves
    call system('gnuplot gplot.txt',ret)
    call system('rm gplot.txt',ret)
    call system('rm xydata.dat',ret)
  end subroutine 

  subroutine gnu_lattice_plot(S,plot_no,title)
    integer, intent(in) :: S(:,:), plot_no
    character(*), intent(in) :: title
    integer :: i, j, ret
    character(30) :: rowfmt, filename

    write(rowfmt, '(A,I4,A)') '(',L,'(1X,I3))' 
    write(filename,'(A,I1,A)') 'set output "plot',plot_no,'.png"'
    
    open(10,access = 'sequential',file = 'Sdata.dat')
    do i = 1,L
      write(10,rowfmt) (S(i,j), j=1,L) ! write spin configuration to file
    enddo
    close(10,status= 'keep')
    
    open(10,access = 'sequential',file = 'matplot.txt')
    ! set output terminal  
    write(10,*) 'set term pngcairo size 640,640 enhanced font "Verdana,10"'
    write(10,*) filename
    write(10,*) 'set border linewidth 0'
    !write(10,*) 'unset key'
    !write(10,*) 'unset colorbox'
    !write(10,*) 'unset tics'
    write(10,*) 'set lmargin screen 0.1'
    write(10,*) 'set rmargin screen 0.9'
    write(10,*) 'set tmargin screen 0.9'
    write(10,*) 'set bmargin screen 0.1'
    write(10,*) 'set palette maxcolors 2'
    write(10,*) 'set palette defined ( -1 "#0066ff", 1 "#ff3300")'
    write(10,*) 'set cbrange [-1:1]'
    write(10,*) 'set cbtics ("+" 1, "-" -1)'
    write(10,*) 'set title "'//TRIM(title)//'"'

    write(10,*) 'set pm3d map'
    write(10,*) 'splot "Sdata.dat" matrix with image'
    close(10,status = 'keep')

    ! now call gnuplot and plot the matrix
    call system('gnuplot matplot.txt',ret)
    call system('rm matplot.txt',ret)
    call system('rm Sdata.dat',ret)
  end subroutine

end module 
