module input_mod
  use Data_Types
  use finer, only: file_ini
  use IO_mod
  use linspace_mod
  use interp1_mod
  implicit none
  private
  public :: read_input

contains

  subroutine read_input(grid, fini, OUT_Path, write_coarse)
    type(Grid_Type),    intent(inout) :: grid
    type(file_ini),     intent(inout) :: fini
    character(len=dl),  intent(out)   :: OUT_Path
    logical,            intent(out)   :: write_coarse

    logical :: ans
    character(len=dl) :: filename, name
    character(len=4) :: indb, indm, indn
    character(len=:), allocatable :: list(:)
    integer :: b, m, n, i, dim(2), error
    integer, allocatable :: num(:)
    real(R8), allocatable :: array(:)
    real(R8) :: cut(2)

    call fini%load(filename='input.ini')

    call fini%get_sections_list(list=list)
    grid%nblocks = 0
    do i = 1, size(list,1)
      if (index(list(i),'PATCHWORK-Block')>0) then
        if (index(list(i),'Tile')==0) grid%nblocks = grid%nblocks + 1
      end if
    end do
    allocate( grid%blk(1:grid%nblocks))

    call fini%get(section_name='PATCHWORK', option_name='outpath', val=OUT_Path, error=error)
    if (error/=0) OUT_Path = './'

    call fini%get(section_name='PATCHWORK', option_name='coarse-mesh', val=write_coarse, error=error)
    if (error/=0) write_coarse = .false.

    ! Block dimensions
    do b = 1, grid%nblocks
      write(indb,'(I4)') b
      call fini%get(section_name='PATCHWORK-Block'//adjustl(indb), option_name='tiles', val=dim(:), error=error)
      grid%blk(b)%m = dim(1)
      grid%blk(b)%n = dim(2)
      allocate( grid%blk(b)%tile( grid%blk(b)%m, grid%blk(b)%n) )
    end do

    ! Complete each block info
    do b = 1, grid%nblocks
      write(*,*) ' Reading block', b
      write(indb,'(I4)') b

      ! ni per tile column
      if (allocated(num)) deallocate(num)
      allocate( num(grid%blk(b)%m) )
      call fini%get(section_name='PATCHWORK-Block'//adjustl(indb), option_name='ni', val=num, error=error)
      do n = 1, grid%blk(b)%n
        do m = 1,grid%blk(b)%m
          grid%blk(b)%tile(m,n)%ni = num(m)
        end do
      end do

      ! nj per tile row
      if (allocated(num)) deallocate(num)
      allocate( num(grid%blk(b)%n) )
      call fini%get(section_name='PATCHWORK-Block'//adjustl(indb), option_name='nj', val=num, error=error)
      do n = 1, grid%blk(b)%n
        do m = 1,grid%blk(b)%m
          grid%blk(b)%tile(m,n)%nj = num(n)
        end do
      end do

      do n = 1,grid%blk(b)%n
        do m = 1,grid%blk(b)%m
          write(*,*) '  Tile', m, 'x', n
          associate( t => grid%blk(b)%tile(m,n) )

          allocate( t%x(t%ni+1, t%nj+1) )
          allocate( t%y(t%ni+1, t%nj+1) )

          ! Default initialization for flags (fixes uninitialized 'connect')
          t%i_str_smooth = ' no'
          t%j_str_smooth = ' no'
          t%smooth_low   = ' no'
          t%smooth_upp   = ' no'
          t%connect      = ' no'
          t%def_jstr_sm_range = .false.

          write(indm,'(I4)') m
          write(indn,'(I4)') n
          name = 'PATCHWORK-Block'//trim(adjustl(indb))//'-Tile'//trim(adjustl(indm))//trim(adjustl(indn))

          call fini%get(section_name=trim(name), option_name='dx', val=t%dx, error=error)
          if (error/=0) t%dx = 1000._R8
          call fini%get(section_name=trim(name), option_name='dy', val=t%dy, error=error)
          if (error/=0) t%dy = 1000._R8
          call fini%get(section_name=trim(name), option_name='sx', val=t%sx, error=error)
          if (error/=0) t%sx = '=='
          call fini%get(section_name=trim(name), option_name='sy', val=t%sy, error=error)
          if (error/=0) t%sy = '=='

          ! i stretch transitions
          call fini%get(section_name=trim(name), option_name='dx2', val=t%dx_s, error=error)
          if (error==0) t%i_str_smooth = 'yes'
          call fini%get(section_name=trim(name), option_name='sx2', val=t%sx_s, error=error)
          if (error/=0) t%sx_s = t%sx
          call fini%get(section_name=trim(name), option_name='i-str-change-dj', val=dim, error=error)
          if (error/=0) then
            t%j_st_ism  = dim(1)
            t%j_end_ism = dim(2)
          else
            t%j_st_ism  = 1
            t%j_end_ism = t%nj + 1
          end if

          ! j stretch transitions
          call fini%get(section_name=trim(name), option_name='dy2', val=t%dy_s, error=error)
          if (error==0) then
            t%j_str_smooth = 'yes'
            call fini%get(section_name=trim(name), option_name='transition-x-range', &
                          val=t%x_jstr_sm_range(:), error=error)
            if (error==0) then
              t%def_jstr_sm_range = .true.
            else
              t%def_jstr_sm_range = .false.
            end if
          end if
          call fini%get(section_name=trim(name), option_name='sy2', val=t%sy_s, error=error)
          if (error/=0) t%sy_s = t%sy
          call fini%get(section_name=trim(name), option_name='sf',  val=t%s_f,  error=error)

          ! Cubic vs straight lower/upper segments
          call fini%get(section_name=trim(name), option_name='cubic-lower-profile', val=ans, error=error)
          if ((error==0) .and. (ans)) t%smooth_low = 'yes'
          call fini%get(section_name=trim(name), option_name='cubic-upper-profile', val=ans, error=error)
          if ((error==0) .and. (ans)) t%smooth_upp = 'yes'

          ! Connection
          call fini%get(section_name=trim(name), option_name='connection-data', val=t%ind_con(:), error=error)
          if (error==0) t%connect = 'yes'

          ! Lower/upper profiles
          call read_profile(fini, trim(name), 'lower-profile', t%xp_low, t%yp_low, t%np_low)
          call read_profile(fini, trim(name), 'upper-profile', t%xp_upp, t%yp_upp, t%np_upp)

          write(*,*) '    -> done'
          end associate
        end do
      end do
    end do

  end subroutine read_input


  ! Reads one lower/upper profile: first tries filename+cut, then falls back
  ! to an inline array given as x1 y1 x2 y2 ...
  subroutine read_profile(fini, section, option_base, xp, yp, npoints)
    type(file_ini),   intent(inout) :: fini
    character(len=*), intent(in)    :: section, option_base
    real(R8), allocatable, intent(out) :: xp(:), yp(:)
    integer,  intent(out)           :: npoints

    character(len=dl) :: filename
    real(R8), allocatable :: array(:)
    real(R8) :: cut(2)
    integer :: error, i, nvals

    call fini%get(section_name=section, option_name=option_base, val=filename, error=error)
    call importdata(filename=filename, error=error, x=xp, y=yp)
    if (error==0) then
      call fini%get(section_name=section, option_name=option_base//'-cut', val=cut, error=error)
      if (error==0) call trim_data(xp, yp, cut(1), cut(2))
      npoints = size(xp)
    else
      nvals = fini%count_values(section_name=section, option_name=option_base)
      allocate(array(nvals))
      call fini%get(section_name=section, option_name=option_base, val=array, error=error)
      npoints = size(array)/2
      allocate(xp(npoints), yp(npoints))
      do i = 1, npoints
        xp(i) = array(2*i-1)
        yp(i) = array(2*i)
      end do
    end if
  end subroutine read_profile


  subroutine trim_data(x, y, x1, x2)
    real(R8), intent(in)    :: x1, x2
    real(R8), intent(inout) :: x(:), y(:)
    real(R8), allocatable   :: xtemp(:), ytemp(:)

    allocate(xtemp(size(x)), ytemp(size(x)))
    xtemp = linspace(x1, x2, size(x))
    call interp1(x, y, xtemp, ytemp)
    x = xtemp
    y = ytemp
  end subroutine trim_data

end module input_mod
