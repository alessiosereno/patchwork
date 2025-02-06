program patchwork
  use Data_Types
  use finer, only: file_ini
  use IO_mod
  use linspace_mod
  use stretching_mod
  use interp1_mod
  use spline_mod
  use cubic_mod
  implicit none

  logical :: trapez_1, trapez_2
  integer :: b, i, j, m, n, bs, ncell
  integer :: js1, js2, cursor_i, cursor_j
  character(8), parameter :: real_form = '(e20.10)'
  character(49), parameter :: mesh_dim = "(' I= ',i3,', J= ',i3,', K= 1, ZONETYPE=Ordered')"
  character(len=dl) :: inputfile, OUT_Path
  logical :: write_coarse
  real(8), allocatable, dimension(:) :: x1, y1, y2, y_line_1, y_line_2, y_line_3, x_line_1, x_line_2
  real(8) :: weight_1, weight_2
  real(8) :: lenght_loc, lenght_low, lenght_upp
  integer :: conn_nj, conn_ni, loca_nj
  type(Grid_Type) :: grid
  type(Point_Type) :: point_1, point_2
  type(file_ini) :: fini

  ! Read input file
  call read_input()

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Block dimensions
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  write(*,*), ' Allocating blocks.'
  do b = 1, grid%nblocks
   grid%blk(b)%ni = 0
   grid%blk(b)%nj = 0
    n = 1
    do m = 1,grid%blk(b)%m
     grid%blk(b)%ni = grid%blk(b)%ni +grid%blk(b)%tile(m,n)%ni
    end do
    m = 1
    do n = 1,grid%blk(b)%n
     grid%blk(b)%nj = grid%blk(b)%nj +grid%blk(b)%tile(m,n)%nj
    end do
    write(*,*) '  Block ', b,grid%blk(b)%ni,'x',grid%blk(b)%nj
    allocate( grid%blk(b)%x(grid%blk(b)%ni+1,grid%blk(b)%nj+1) )
    allocate( grid%blk(b)%y(grid%blk(b)%ni+1,grid%blk(b)%nj+1) )
  end do
  ! - - - - - - - - - - - - - - - - - - - - - - - - -

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: x-constant lines
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, grid%nblocks
    do n = 1,grid%blk(b)%n
      do m = 1,grid%blk(b)%m

        if ( grid%blk(b)%tile(m,n)%i_str_smooth /= 'yes' ) then

          do j = 1,grid%blk(b)%tile(m,n)%nj+1

           grid%blk(b)%tile(m,n)%x(:,j) = &
            stretch ( start = grid%blk(b)%tile(m,n)%xp_low(1),   &
                      end   = grid%blk(b)%tile(m,n)%xp_low(      &
                             grid%blk(b)%tile(m,n)%np_low),     &
                      n     = grid%blk(b)%tile(m,n)%ni+1,        &
                      dir   = grid%blk(b)%tile(m,n)%sx,          &
                      delta = grid%blk(b)%tile(m,n)%dx ) 
          end do

        elseif ( grid%blk(b)%tile(m,n)%i_str_smooth == 'yes' ) then
            
          allocate( x_line_1( grid%blk(b)%tile(m,n)%ni+1 ) )
          allocate( x_line_2( grid%blk(b)%tile(m,n)%ni+1 ) )
          
          do j = 1,grid%blk(b)%tile(m,n)%nj+1

            x_line_1 = stretch ( start = grid%blk(b)%tile(m,n)%xp_low(1),  &
                                 end   = grid%blk(b)%tile(m,n)%xp_low(     &
                                         grid%blk(b)%tile(m,n)%np_low),    &
                                 n     = grid%blk(b)%tile(m,n)%ni+1,       &
                                 dir   = grid%blk(b)%tile(m,n)%sx,         &
                                 delta = grid%blk(b)%tile(m,n)%dx ) 

            x_line_2 = stretch ( start = grid%blk(b)%tile(m,n)%xp_low(1),  &
                                 end   = grid%blk(b)%tile(m,n)%xp_low(     &
                                         grid%blk(b)%tile(m,n)%np_low),    &
                                 n     = grid%blk(b)%tile(m,n)%ni+1,       &
                                 dir   = grid%blk(b)%tile(m,n)%sx_s,       &
                                 delta = grid%blk(b)%tile(m,n)%dx_s ) 
            
            weight_1 = float ( ( j -grid%blk(b)%tile(m,n)%j_st_ism ) ) / &
                       float ( grid%blk(b)%tile(m,n)%j_end_ism -grid%blk(b)%tile(m,n)%j_st_ism ) 
            weight_1 = 2d0*weight_1**3 - 3d0*weight_1**2 + 1d0
            if ( j <grid%blk(b)%tile(m,n)%j_st_ism  ) weight_1 = 1d0
            if ( j >grid%blk(b)%tile(m,n)%j_end_ism ) weight_1 = 0d0
            weight_2 = 1d0 - weight_1
           grid%blk(b)%tile(m,n)%x(:,j) = x_line_1 * weight_1 + x_line_2 * weight_2

          end do

          deallocate( x_line_1, x_line_2 )

        end if
      
      end do
    end do
  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: y discretization
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, grid%nblocks
   
    ! Overall block dim, to be updated after each tile
   grid%blk(b)%ni = 0
   grid%blk(b)%nj = 0
    
    do n = 1,grid%blk(b)%n
      do m = 1,grid%blk(b)%m

        allocate( x1( grid%blk(b)%tile(m,n)%ni+1 ) )
        allocate( y1( grid%blk(b)%tile(m,n)%ni+1 ) )
        allocate( y2( grid%blk(b)%tile(m,n)%ni+1 ) )
        
        ! interpolate lower and upper boundaries of the tile
        x1 = grid%blk(b)%tile(m,n)%x(:,1)
        if ( grid%blk(b)%tile(m,n)%np_low < 10 ) then
          call interp1( grid%blk(b)%tile(m,n)%xp_low,grid%blk(b)%tile(m,n)%yp_low, x1, y1 )
        else
          y1 = spline( grid%blk(b)%tile(m,n)%xp_low,grid%blk(b)%tile(m,n)%yp_low, x1 )
        end if
        x1 = grid%blk(b)%tile(m,n)%x(:,grid%blk(b)%tile(m,n)%nj+1)
        if ( grid%blk(b)%tile(m,n)%np_upp < 10 ) then
          call interp1( grid%blk(b)%tile(m,n)%xp_upp,grid%blk(b)%tile(m,n)%yp_upp, x1, y2 )
        else
          y2 = spline( grid%blk(b)%tile(m,n)%xp_upp,grid%blk(b)%tile(m,n)%yp_upp, x1 )
        end if

        ! option to add a cubic transition from first to last point
        if ( grid%blk(b)%tile(m,n)%smooth_low == 'yes' ) then
          point_1%x = x1(1)
          point_2%x = x1(size(x1))
          point_1%y = grid%blk(b)%tile(m,n)%yp_low(1)
          point_2%y = grid%blk(b)%tile(m,n)%yp_low(size(grid%blk(b)%tile(m,n)%yp_low))
          call cubic( point_1%x, point_1%y, point_2%x, point_2%y, x1, y1 )
        end if

        if ( grid%blk(b)%tile(m,n)%smooth_upp == 'yes' ) then
          point_1%x = x1(1)
          point_2%x = x1(size(x1))
          point_1%y = grid%blk(b)%tile(m,n)%yp_upp(1)
          point_2%y = grid%blk(b)%tile(m,n)%yp_upp(size(grid%blk(b)%tile(m,n)%yp_upp))
          call cubic( point_1%x, point_1%y, point_2%x, point_2%y, x1, y2 )
        end if

        if ( grid%blk(b)%tile(m,n)%j_str_smooth /= 'yes' .and. &
            grid%blk(b)%tile(m,n)%connect /= 'yes' ) then

          do i = 1,grid%blk(b)%tile(m,n)%ni + 1

           grid%blk(b)%tile(m,n)%y(i,:) =                 &
            stretch ( start = y1(i),                  &
                      end   = y2(i),                  &
                      n     = grid%blk(b)%tile(m,n)%nj+1,  &
                      dir   = grid%blk(b)%tile(m,n)%sy,    &
                      delta = grid%blk(b)%tile(m,n)%dy )
          end do
    
        elseif ( grid%blk(b)%tile(m,n)%j_str_smooth == 'yes' .and. &
                grid%blk(b)%tile(m,n)%connect /= 'yes' ) then

          allocate( y_line_1( grid%blk(b)%tile(m,n)%nj+1 ) )
          allocate( y_line_2( grid%blk(b)%tile(m,n)%nj+1 ) )

          do i = 1,grid%blk(b)%tile(m,n)%ni + 1
            y_line_1 = stretch ( start = y1(i),                  &
                                 end   = y2(i),                  &
                                 n     = grid%blk(b)%tile(m,n)%nj+1,  &
                                 dir   = grid%blk(b)%tile(m,n)%sy,    &
                                 delta = grid%blk(b)%tile(m,n)%dy ) 

            y_line_2 = stretch ( start = y1(i),                  &
                                 end   = y2(i),                  &
                                 n     = grid%blk(b)%tile(m,n)%nj+1,  &
                                 dir   = grid%blk(b)%tile(m,n)%sy_s,  &
                                 delta = grid%blk(b)%tile(m,n)%dy_s ) 

            weight_1 = ( x1(i) - x1(1) ) / ( x1(size(x1)) - x1(1) )
            if ( weight_1<0.1 ) then
              grid%blk(b)%tile(m,n)%y(i,:) = y_line_1
            elseif ( weight_1>0.9 ) then
              grid%blk(b)%tile(m,n)%y(i,:) = y_line_2
            else
              weight_1 = (2d0*weight_1**3 - 3d0*weight_1**2 + 0.54*weight_1)/0.512d0 + 0.949219d0
              weight_2 = 1d0 - weight_1
              grid%blk(b)%tile(m,n)%y(i,:) = y_line_1 * weight_1 + y_line_2 * weight_2
            endif

          end do

          deallocate( y_line_1, y_line_2 )
    
        elseif ( grid%blk(b)%tile(m,n)%connect=='yes' ) then

          ! work arrays
          allocate( y_line_1( grid%blk(b)%tile(m,n)%nj+1 ) )
          allocate( y_line_2( grid%blk(b)%tile(m,n)%nj+1 ) )
          allocate( y_line_3( grid%blk(b)%tile(m,n)%nj+1 ) )

          ! Connection info
          bs = grid%blk(b)%tile(m,n)%ind_con(1)
          js1 = grid%blk(b)%tile(m,n)%ind_con(2)
          js2 = grid%blk(b)%tile(m,n)%ind_con(3)
          conn_nj = grid%blk(bs)%nj
          conn_ni = grid%blk(bs)%ni
          loca_nj = grid%blk(b)%tile(m,n)%nj

          do i = 1,grid%blk(b)%tile(m,n)%ni + 1
            
            ! Discretized line from connected block
            y_line_1(1:loca_nj+1) = grid%blk(bs)%y(conn_ni+1,js1:js2)
            !y_line_3(1:conn_nj+1) = linspace(y_line_1(1), y_line_1(conn_nj+1), conn_nj+1)
            weight_1 = ( x1(i) - x1(1) ) / ( x1(size(x1)) - x1(1) )
            weight_1 = 2d0*weight_1**3 - 3d0*weight_1**2 + 1d0
            weight_1 = max( weight_1, 0d0 )
            weight_2 = 1d0 - weight_1
            !y_line_1(1:conn_nj+1) = y_line_1(1:conn_nj+1)*weight_1 + &
            !                        y_line_3(1:conn_nj+1)*weight_2
          
            ! Reduce/enlarge y_line_1 to take into account y_upp of local wall
            y_line_1(1:loca_nj+1) = ( y_line_1(1:loca_nj+1) - y_line_1(1) ) &
                                  * ( y2(i) - y1(i) )/( y2(1) - y1(1) )  &
                                  + y_line_1(1)

            ! Right y-line for this tile
            y_line_2 = stretch ( start = y1(i),                 &
                                 end   = y2(i),                 &
                                 n     = grid%blk(b)%tile(m,n)%nj+1, &
                                 dir   = grid%blk(b)%tile(m,n)%sy,   &
                                 delta = grid%blk(b)%tile(m,n)%dy    ) 

            ! Smoothing from left to right discretization
            grid%blk(b)%tile(m,n)%y(i,:) = y_line_1 * weight_1 + y_line_2 * weight_2

          end do ! tile i-loop

          deallocate( y_line_1, y_line_2, y_line_3 )

        end if ! standard / smoothing / connection if
    
        deallocate( x1, y1, y2 )

        ! Adjust x-discretization for trapezoidal-shaped tiles
        trapez_1 = grid%blk(b)%tile(m,n)%xp_low(1) /= grid%blk(b)%tile(m,n)%xp_upp(1)
        trapez_2 = grid%blk(b)%tile(m,n)%xp_low(grid%blk(b)%tile(m,n)%np_low) /= &
                   grid%blk(b)%tile(m,n)%xp_upp(grid%blk(b)%tile(m,n)%np_upp)

        if ( trapez_1 .or. trapez_2 ) then

          lenght_low = grid%blk(b)%tile(m,n)%xp_low(grid%blk(b)%tile(m,n)%np_low) &
                      -grid%blk(b)%tile(m,n)%xp_low(1)
          lenght_upp = grid%blk(b)%tile(m,n)%xp_upp(grid%blk(b)%tile(m,n)%np_upp) &
                      -grid%blk(b)%tile(m,n)%xp_upp(1)

          do j = 1,grid%blk(b)%tile(m,n)%nj + 1

            weight_1 = grid%blk(b)%tile(m,n)%y(1,grid%blk(b)%tile(m,n)%nj+1) -grid%blk(b)%tile(m,n)%y(1,1)
            weight_1 = ( grid%blk(b)%tile(m,n)%y(1,j) -grid%blk(b)%tile(m,n)%y(1,1) ) / weight_1
            weight_2 = 1d0 - weight_1
            lenght_loc = weight_2 * lenght_low + weight_1 * lenght_upp   
            
            if ( trapez_1 .and. .not. trapez_2 ) then
             grid%blk(b)%tile(m,n)%x(:,j) = ( grid%blk(b)%tile(m,n)%x(:,j) -grid%blk(b)%tile(m,n)%x(1,j) ) &
                                      / ( grid%blk(b)%tile(m,n)%x(grid%blk(b)%tile(m,n)%ni+1,j) &
                                      -  grid%blk(b)%tile(m,n)%x(1,j) ) &
                                      * lenght_loc +grid%blk(b)%tile(m,n)%x(1,j) &
                                      + (lenght_low - lenght_upp) * weight_1

            elseif ( trapez_2 .and. .not. trapez_1 ) then  
             grid%blk(b)%tile(m,n)%x(:,j) = ( grid%blk(b)%tile(m,n)%x(:,j) -grid%blk(b)%tile(m,n)%x(1,j) ) &
                                      / ( grid%blk(b)%tile(m,n)%x(grid%blk(b)%tile(m,n)%ni+1,j) &
                                      -  grid%blk(b)%tile(m,n)%x(1,j) ) &
                                      * lenght_loc +grid%blk(b)%tile(m,n)%x(1,j)

            elseif ( trapez_1 .and. trapez_2 ) then
             grid%blk(b)%tile(m,n)%x(:,j) = ( grid%blk(b)%tile(m,n)%x(:,j) -grid%blk(b)%tile(m,n)%x(1,j) ) &
                                      / ( grid%blk(b)%tile(m,n)%x(grid%blk(b)%tile(m,n)%ni+1,j) &
                                      -  grid%blk(b)%tile(m,n)%x(1,j) ) &
                                      * lenght_loc +grid%blk(b)%tile(m,n)%x(1,j) &
                                      + (lenght_low - lenght_upp) * weight_1 /2d0

            end if
        
          end do
        
        end if

        if ( n == 1 ) then
          grid%blk(b)%ni = grid%blk(b)%ni +grid%blk(b)%tile(m,n)%ni
          write(*,*) ' Updating i-dimension of block', b, '->',grid%blk(b)%ni
        end if

      end do ! tile m loop

      grid%blk(b)%nj = grid%blk(b)%nj +grid%blk(b)%tile(1,n)%nj
      write(*,*) ' Updating j-dimension of block', b, '->',grid%blk(b)%nj
    
    end do ! tile n loop

    write(*,*) ' Assembling block', b
    cursor_j = 1
    do n = 1,grid%blk(b)%n

      cursor_i = 1
      do m = 1,grid%blk(b)%m

        grid%blk(b)%x( cursor_i : cursor_i +grid%blk(b)%tile(m,n)%ni,   &
                       cursor_j : cursor_j +grid%blk(b)%tile(m,n)%nj) = &
                       grid%blk(b)%tile(m,n)%x

        grid%blk(b)%y( cursor_i : cursor_i +grid%blk(b)%tile(m,n)%ni,   &
                       cursor_j : cursor_j +grid%blk(b)%tile(m,n)%nj) = &
                       grid%blk(b)%tile(m,n)%y

        cursor_i = cursor_i +grid%blk(b)%tile(m,n)%ni

      end do
      cursor_j = cursor_j + grid%blk(b)%tile(m-1,n)%nj

    end do
    write(*,*) '    -> done'

  end do ! block loop

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Writing the mesh
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  open(10,file=trim(OUT_Path)//'mesh.dat',status='unknown')
  write(10,*) 'TITLE     = "Mesh"'
  write(10,*) 'VARIABLES = "x"'
  write(10,*) '"y"'

  do b = 1, grid%nblocks

    write(10,*)'ZONE T="Block ', b,'"'
    
    write(10,mesh_dim)grid%blk(b)%ni+1,grid%blk(b)%nj+1
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'DT=(SINGLE SINGLE)'

    do j = 1,grid%blk(b)%nj + 1
      do i = 1,grid%blk(b)%ni + 1
        write (10,real_form) grid%blk(b)%x(i,j)
      end do
    end do

    do j = 1,grid%blk(b)%nj + 1
      do i = 1,grid%blk(b)%ni +1
        write (10,real_form) grid%blk(b)%y(i,j)
      end do
    end do

  end do
  close(10)

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Total number of cells
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  ncell = 0
  do b = 1, grid%nblocks
    ncell = ncell + grid%blk(b)%ni*grid%blk(b)%nj
  end do
  write(*,*) ' Total number of cells', ncell

  if (write_coarse) then
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Writing the coarse meshes
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
    write(*,*) ' 2nd level mesh with: ', ncell/4, 'cells'
    open(10,file=trim(OUT_Path)//'mesh2.dat',status='unknown')
    write(10,*) 'TITLE     = "Mesh"'
    write(10,*) 'VARIABLES = "x"'
    write(10,*) '"y"'

    do b = 1, grid%nblocks

      write(10,*)'ZONE T="Block ', b,'"'
      
      write(10,mesh_dim)grid%blk(b)%ni/2+1,grid%blk(b)%nj/2+1
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'DT=(SINGLE SINGLE)'

      do j = 1,grid%blk(b)%nj + 1, 2
        do i = 1,grid%blk(b)%ni + 1, 2
          write (10,real_form) grid%blk(b)%x(i,j)
        end do
      end do

      do j = 1,grid%blk(b)%nj + 1, 2
        do i = 1,grid%blk(b)%ni + 1, 2
          write (10,real_form) grid%blk(b)%y(i,j)
        end do
      end do

    end do
    close(10)

    write(*,*) ' 3nd level mesh with: ', ncell/16, 'cells'
    open(10,file=trim(OUT_Path)//'mesh3.dat',status='unknown')
    write(10,*) 'TITLE     = "Mesh"'
    write(10,*) 'VARIABLES = "x"'
    write(10,*) '"y"'

    do b = 1, grid%nblocks

      write(10,*)'ZONE T="Block ', b,'"'
      
      write(10,mesh_dim)grid%blk(b)%ni/4+1,grid%blk(b)%nj/4+1
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'DT=(SINGLE SINGLE)'

      do j = 1,grid%blk(b)%nj + 1, 4
        do i = 1,grid%blk(b)%ni + 1, 4
          write (10,real_form) grid%blk(b)%x(i,j)
        end do
      end do

      do j = 1,grid%blk(b)%nj + 1, 4
        do i = 1,grid%blk(b)%ni + 1, 4
          write (10,real_form) grid%blk(b)%y(i,j)
        end do
      end do

    end do
    close(10)
  end if

contains

  subroutine read_input()
    implicit none
    logical :: ans
    character(len=dl) :: includefile, filename, name
    character(len=4) :: indb, indm, indn
    character(len=:), allocatable :: list(:) !< List of all section names.
    integer :: b, dim(2), error
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
    enddo
    allocate( grid%blk(1:grid%nblocks))

    ! Assignement of the OUT_Path
    call fini%get(section_name='PATCHWORK', option_name='outpath', val=OUT_Path, error=error)
    if (error/=0) OUT_Path = './'

    ! Write coarse meshes
    call fini%get(section_name='PATCHWORK', option_name='coarse-mesh', val=write_coarse, error=error)
    if (error/=0) write_coarse = .false.

    ! Blocks general info
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

      ! ni of tiles ('m' values)
      if (allocated(num)) deallocate(num)
      allocate( num(grid%blk(b)%m) )
      call fini%get(section_name='PATCHWORK-Block'//adjustl(indb), option_name='ni', val=num, error=error)
      do n = 1, grid%blk(b)%n
        do m = 1,grid%blk(b)%m
          grid%blk(b)%tile(m,n)%ni = num(m)
        end do
      end do

      ! nj of tiles ('n' values)
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
          ! tile allocation
          allocate( grid%blk(b)%tile(m,n)%x(grid%blk(b)%tile(m,n)%ni+1,grid%blk(b)%tile(m,n)%nj+1) )
          allocate( grid%blk(b)%tile(m,n)%y(grid%blk(b)%tile(m,n)%ni+1,grid%blk(b)%tile(m,n)%nj+1) )

          write(indm,'(I4)') m
          write(indn,'(I4)') n
          name = 'PATCHWORK-Block'//trim(adjustl(indb))//'-Tile'//trim(adjustl(indm))//trim(adjustl(indn))
          call fini%get(section_name=trim(name), &
                        option_name='dx', val=grid%blk(b)%tile(m,n)%dx, error=error)
          if (error/=0) grid%blk(b)%tile(m,n)%dx = 1000.
          call fini%get(section_name=trim(name), &
                        option_name='dy', val=grid%blk(b)%tile(m,n)%dy, error=error)
          if (error/=0) grid%blk(b)%tile(m,n)%dy = 1000.
          call fini%get(section_name=trim(name), &
                        option_name='sx', val=grid%blk(b)%tile(m,n)%sx, error=error)
          if (error/=0) grid%blk(b)%tile(m,n)%sx = '=='
          call fini%get(section_name=trim(name), &
                        option_name='sy', val=grid%blk(b)%tile(m,n)%sy, error=error)
          if (error/=0) grid%blk(b)%tile(m,n)%sy = '=='

          ! i stretch transitions
          grid%blk(b)%tile(m,n)%i_str_smooth = ' no' ! default
          call fini%get(section_name=trim(name), &
                        option_name='dx2', val=grid%blk(b)%tile(m,n)%dx_s, error=error)
          if (error==0) grid%blk(b)%tile(m,n)%i_str_smooth = 'yes'
          call fini%get(section_name=trim(name), &
                        option_name='sx2', val=grid%blk(b)%tile(m,n)%sx_s, error=error)
          if (error/=0) grid%blk(b)%tile(m,n)%sx_s = grid%blk(b)%tile(m,n)%sx
          call fini%get(section_name=trim(name), &
                        option_name='i-str-change-dj', val=dim, error=error)
          if (error/=0) then
            grid%blk(b)%tile(m,n)%j_st_ism = dim(1)
            grid%blk(b)%tile(m,n)%j_end_ism = dim(2)
          else
            grid%blk(b)%tile(m,n)%j_st_ism = 1
            grid%blk(b)%tile(m,n)%j_end_ism = grid%blk(b)%tile(m,n)%nj + 1
          end if

          ! j stretch transitions
          grid%blk(b)%tile(m,n)%j_str_smooth = ' no' ! default
          call fini%get(section_name=trim(name), &
                        option_name='dy2', val=grid%blk(b)%tile(m,n)%dy_s, error=error)
          if (error==0) grid%blk(b)%tile(m,n)%j_str_smooth = 'yes'
          call fini%get(section_name=trim(name), &
                        option_name='sy2', val=grid%blk(b)%tile(m,n)%sy_s, error=error)
          if (error/=0) grid%blk(b)%tile(m,n)%sy_s = grid%blk(b)%tile(m,n)%sy
          call fini%get(section_name=trim(name), &
                        option_name='sf', val=grid%blk(b)%tile(m,n)%s_f, error=error)

          ! Cubic vs straight lower/upper segments
          grid%blk(b)%tile(m,n)%smooth_low = ' no' ! default
          grid%blk(b)%tile(m,n)%smooth_upp = ' no' ! default
          call fini%get(section_name=trim(name), &
                        option_name='cubic-lower-segment', val=ans, error=error)
          if ((error==0) .and. (ans)) grid%blk(b)%tile(m,n)%smooth_low = 'yes'
          call fini%get(section_name=trim(name), &
                        option_name='cubic-upper-segment', val=ans, error=error)
          if ((error==0) .and. (ans)) grid%blk(b)%tile(m,n)%smooth_upp = 'yes'

          ! Connection
          call fini%get(section_name=trim(name), &
                        option_name='connection-data', val=grid%blk(b)%tile(m,n)%ind_con(:), error=error)
          if (error==0) grid%blk(b)%tile(m,n)%connect = 'yes'

          ! Lower segment
          call fini%get(section_name=trim(name), &
                        option_name='lower-profile', val=filename, error=error)
          call importdata(filename=filename, error=error, &
                          x=grid%blk(b)%tile(m,n)%xp_low, &
                          y=grid%blk(b)%tile(m,n)%yp_low)
          if (error==0) then
            call fini%get(section_name=trim(name), &
                          option_name='lower-profile-cut', val=cut, error=error)
            if (error==0) then
              call trim_data(x=grid%blk(b)%tile(m,n)%xp_low, &
                             y=grid%blk(b)%tile(m,n)%yp_low,x1=cut(1),x2=cut(2))
            end if
            grid%blk(b)%tile(m,n)%np_low = size(grid%blk(b)%tile(m,n)%xp_low)
          else
            ! Array explicitly given
            if (allocated(array)) deallocate(array)
            allocate(array(1:fini%count_values(section_name=trim(name), option_name='lower-profile')))
            call fini%get(section_name=trim(name), option_name='lower-profile', val=array, error=error)
            grid%blk(b)%tile(m,n)%np_low = size(array)/2
            allocate(grid%blk(b)%tile(m,n)%xp_low(size(array)/2))
            allocate(grid%blk(b)%tile(m,n)%yp_low(size(array)/2))
            do i = 1, size(array)/2
              grid%blk(b)%tile(m,n)%xp_low(i) = array(2*i-1)
              grid%blk(b)%tile(m,n)%yp_low(i) = array(2*i)
            end do
          endif

          ! Upper segment
          call fini%get(section_name=trim(name), &
                        option_name='upper-profile', val=filename, error=error)
          call importdata(filename=filename, error=error, &
                          x=grid%blk(b)%tile(m,n)%xp_upp, &
                          y=grid%blk(b)%tile(m,n)%yp_upp)
          if (error==0) then
            call fini%get(section_name=trim(name), &
                          option_name='upper-profile-cut', val=cut, error=error)
            if (error==0) then
              call trim_data(x=grid%blk(b)%tile(m,n)%xp_upp, &
                             y=grid%blk(b)%tile(m,n)%yp_upp,x1=cut(1),x2=cut(2))
            end if
            grid%blk(b)%tile(m,n)%np_upp = size(grid%blk(b)%tile(m,n)%xp_upp)
          else
            ! Array explicitly given
            if (allocated(array)) deallocate(array)
            allocate(array(1:fini%count_values(section_name=trim(name), option_name='upper-profile')))
            call fini%get(section_name=trim(name), option_name='upper-profile', val=array, error=error)
            grid%blk(b)%tile(m,n)%np_upp = size(array)/2
            allocate(grid%blk(b)%tile(m,n)%xp_upp(size(array)/2))
            allocate(grid%blk(b)%tile(m,n)%yp_upp(size(array)/2))
            do i = 1, size(array)/2
              grid%blk(b)%tile(m,n)%xp_upp(i) = array(2*i-1)
              grid%blk(b)%tile(m,n)%yp_upp(i) = array(2*i)
            end do
          endif

          write(*,*) '    -> done'
        end do
      end do ! tile processing
    end do ! block processing

  end subroutine read_input

  subroutine trim_data(x,y,x1,x2)
    implicit none
    real(R8), intent(in) :: x1, x2
    real(R8), intent(inout) :: x(:), y(:)
    real(R8), allocatable :: xtemp(:), ytemp(:)
    integer(I4) :: i

    allocate(xtemp(size(x)))
    allocate(ytemp(size(x)))
    xtemp = linspace(x1, x2, size(x))
    call interp1(x, y, xtemp, ytemp)
    x = xtemp
    y = ytemp

  end subroutine trim_data

end program patchwork
