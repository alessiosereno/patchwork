module IO_mod
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  public :: extract_path

  integer(kind=I4), parameter, public :: dl=500

contains

  subroutine importdata(filename,x,y,error)
    character(len=dl),intent(in) :: filename
    integer(I4), intent(out) :: error
    integer(I4) :: fid    
    real(r8),allocatable,intent(inout) :: x(:), y(:)
    real(r8),dimension(100000) :: xtemp,ytemp 
    integer(I4) :: i

    open(newunit=fid,file=trim(filename),status='old',iostat=error)
    if (error/=0) return
    do i = 1, 100000
      read(fid,*,end=100) xtemp(i),ytemp(i)
    end do
    100 continue

    allocate(x(i-1),y(i-1))

    x = xtemp(1:i-1)
    y = ytemp(1:i-1)

    close(fid)
  end subroutine importdata


  subroutine readarray(filename,x)
    character(len=dl),intent(in)          :: filename
    integer(I4)                       :: fid    
    real(r8),allocatable,intent(inout)    :: x(:)
    real(r8),dimension(100000)            :: xtemp
    integer(I4)                       :: i

    open(newunit=fid,file=trim(filename),status='old')
    do i = 1, 100000
      read(fid,*,end=100) xtemp(i)
    end do
    100  continue

    allocate(x(i-1))

    x = xtemp(1:i-1)

    close(fid)
  end subroutine readarray

  subroutine writearray(filename,x)
    character(len=dl),intent(in)          :: filename
    integer(I4)                       :: fid    
    real(r8),dimension(:),intent(in)      :: x(:)
    integer(I4)                       :: i, n

    n = size(x)

    open(newunit=fid,file=filename,status='unknown')

    do i = 1, n
      write(fid,*) x(i)
    end do

    close(fid)
  end subroutine writearray


  subroutine fliparray(x)
    real(r8),dimension(:),intent(inout)   :: x(:)
    integer(I4)                       :: i, n
    real(r8),allocatable                  :: xflip(:)

    n = size(x)
    allocate(xflip(n))

    do i = 1, n
      xflip(i) = x(n-i+1)
    end do

    x = xflip
  
  end subroutine fliparray


  pure function extract_path(inpath) result(outpath)
    implicit none
    integer :: index_start
    character(len=dl), intent(in)    :: inpath
    character(len=dl)                :: outpath

    ! Find the last occurrence of '/' in the string
    index_start = len_trim(inpath)
    do while (index_start > 0)
      if (inpath(index_start:index_start) == '/') exit
      index_start = index_start - 1
    end do
    ! Extract the folder path
    outpath = inpath(1:index_start)
    if (outpath=='') outpath = './'

  end function extract_path

 END MODULE io_mod

