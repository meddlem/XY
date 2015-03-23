module plotroutines
  use constants
  implicit none
  private
  public :: line_plot, write_lattice, close_lattice_plot, animate_lattice

contains
  subroutine animate_lattice(title)
    character(*), intent(in) :: title
    integer :: ret
    
    ! creates fifo pipe: plotfifo.dat
    call system("rm -f plotfifo.dat; mkfifo plotfifo.dat",ret)     
    
    ! create a gnuplot command file
    open(10,access = 'sequential',file = 'matplot.plt')
      write(10,*) 'set term wxt' !x11 term for better performance
      write(10,*) 'set border linewidth 0'
      write(10,*) 'set lmargin screen 0.1'
      write(10,*) 'set rmargin screen 0.9'
      write(10,*) 'set tmargin screen 0.9'
      write(10,*) 'set bmargin screen 0.1'
      write(10,*) 'unset tics'
      write(10,*) 'unset key'
      write(10,'(A,I3,A)') 'set xrange [0:', L+1, ']'
      write(10,'(A,I3,A)') 'set yrange [0:', L+1, ']'
      write(10,*) 'set title "'//TRIM(title)//'"'
      write(10,*) 'set pm3d map'
      write(10,*) 'load "loop.plt"'
    close(10)
    
    ! create plot/animate instruction
    open(10,access = 'sequential', file = 'loop.plt')
      write(10,*) 'plot "< cat plotfifo.dat" \'
      write(10,*) 'with vectors head size 0.1,20,60 filled'
      write(10,*) 'pause 0.3'
      write(10,*) 'reread'
    close(10)
    
    ! now fork instance of gnuplot to plot/animate the lattice
    call system("gnuplot matplot.plt &",ret)
  end subroutine
  
  subroutine write_lattice(S)
    real(dp), intent(in) :: S(:,:,:)
    integer :: i, j
    character(30) :: rowfmt
    write(rowfmt, '(A)') '(I3,1X,I3,2X,F6.3,2X,F6.3)' 
    
    ! write to pipe 
    open(11,access = 'sequential',status = 'replace',file = 'plotfifo.dat')
      do i = 1,L
        do j = 1,L
          write(11,rowfmt) i, j, 0.45_dp*S(1,i,j), 0.45_dp*S(2,i,j) 
        enddo
      enddo
    close(11)
  end subroutine

  subroutine close_lattice_plot()
    call system('pkill gnuplot')
    call system('rm -f plotfifo.dat')
  end subroutine

  subroutine line_plot(x,y1,xlabel,ylabel,label1,title,plot_no,y2,label2)
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

    xrange = [xmin-(xmax-xmin)*0.1_dp, xmax+(xmax-xmin)*0.1_dp]
    yrange = [ymin-(ymax-ymin)*0.1_dp, ymax+(ymax-ymin)*0.1_dp]

    open(10,access = 'sequential',file = 'xydata.dat')
      do i=1,m
        if (present(y2)) then
          write(10,*) x(i),y1(i),y2(i) ! write datapoints to file
        else
          write(10,*) x(i),y1(i) ! write datapoints to file
        endif
      enddo
    close(10)
    
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
      
      if (m > 0) then
        write(10,*) 'plot "xydata.dat" using 1:2 with line ls 10 t "", \'
        write(10,*) &
        ' "xydata.dat" using 1:2 with line ls 1 t "'//TRIM(label1)//'", \'
        if (present(y2)) then
          write(10,*) &
          ' "xydata.dat" using 1:3 with line ls 2 t "'//TRIM(label2)//'"'
        endif
      endif
    close(10)

    ! now call gnuplot and plot the curves
    call system('gnuplot gplot.txt',ret)
    call system('rm gplot.txt; rm xydata.dat',ret)
  end subroutine 
end module 
