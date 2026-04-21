module IO_mod
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

  integer(kind=I4), parameter, public :: dl = 500

contains

  subroutine importdata(filename, x, y, error)
    character(len=dl),    intent(in)  :: filename
    integer(I4),          intent(out) :: error
    real(R8), allocatable, intent(inout) :: x(:), y(:)
    real(R8), dimension(100000) :: xtemp, ytemp
    integer(I4) :: fid, i

    open(newunit=fid, file=trim(filename), status='old', iostat=error)
    if (error/=0) return
    do i = 1, 100000
      read(fid, *, end=100) xtemp(i), ytemp(i)
    end do
100 continue

    allocate(x(i-1), y(i-1))
    x = xtemp(1:i-1)
    y = ytemp(1:i-1)

    close(fid)
  end subroutine importdata

end module IO_mod
