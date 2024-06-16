module cubic_mod
  implicit none
  contains
  
  subroutine cubic ( x1, y1, x2, y2, x_eval, y_eval )
    use linspace_mod
    implicit none
    real(8), intent(in) :: x1, y1, x2, y2
    real(8), intent(in), dimension(:) :: x_eval
    real(8), intent(inout), dimension(:) :: y_eval
    ! Local
    real(8), dimension(:), allocatable :: x

    allocate( x(size(x_eval)) )
    x = ( x_eval - x_eval(1) ) / ( x_eval(size(x_eval)) - x_eval(1) )
    if ( y1 > y2 ) then
      print*, 'check'
      y_eval = 2d0 * x**3 - 3d0* x**2 + 1d0
    else
      y_eval = -2d0 * x**3 + 3d0* x**2
    end if
    y_eval = y_eval * abs(y2 - y1) 
    y_eval = y_eval - (y_eval(1)-y1)
    deallocate(x)
  end subroutine cubic
end module cubic_mod