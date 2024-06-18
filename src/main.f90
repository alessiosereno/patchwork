program patchwork
  use Data_Types
  use linspace_mod
  use stretching_mod
  use interp1_mod
  use cubic_mod
  implicit none

  logical :: trapez_1, trapez_2
  integer :: b, i, j, h, m, n, nblocks, bs, ms, ns, ncell
  integer :: cursor_i, cursor_j, file_line
  character(8), parameter :: real_form = '(e20.10)'
  character(49), parameter :: mesh_dim = "(' I= ',i3,', J= ',i3,', K= 1, ZONETYPE=Ordered')"
  real(8), allocatable, dimension(:) :: x1, y1, y2, y_line_1, y_line_2, y_line_3, x_line_1, x_line_2
  real(8) :: weight_1, weight_2
  real(8) :: lenght_loc, lenght_low, lenght_upp
  integer :: conn_nj, conn_ni, loca_nj
  type(Block_Type), allocatable, dimension(:) :: blk
  type(Point_Type) :: point_1, point_2

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Read input file
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  open(10,file='patch.inp')
    ! Blocks general info
    read(10,*) nblocks
    allocate( blk(nblocks) )
    do b = 1, nblocks
      read(10,*) blk(b)%m, blk(b)%n
      allocate( blk(b)%tile(blk(b)%m, blk(b)%n) )
    end do

    ! Read block tiles
    do b = 1, nblocks
      ! Header
      read(10,*) ; read(10,*) ; read(10,*) ; read(10,*)
      write(*,*) ' Reading block', b
      do n = 1, blk(b)%n
        do m = 1, blk(b)%m
          write(*,*) '  Tile', m, 'x', n
          read(10,*) ; read(10,*)
          read(10,*) blk(b)%tile(m,n)%ni
          read(10,*) blk(b)%tile(m,n)%nj
          allocate( blk(b)%tile(m,n)%x(blk(b)%tile(m,n)%ni+1,blk(b)%tile(m,n)%nj+1) )
          allocate( blk(b)%tile(m,n)%y(blk(b)%tile(m,n)%ni+1,blk(b)%tile(m,n)%nj+1) )
          read(10,*) blk(b)%tile(m,n)%dx
          read(10,*) blk(b)%tile(m,n)%dy
          read(10,*) blk(b)%tile(m,n)%sx
          read(10,*) blk(b)%tile(m,n)%sy

          read(10,*) blk(b)%tile(m,n)%i_str_smooth
          if ( blk(b)%tile(m,n)%i_str_smooth == 'yes' ) then
            read(10,*) blk(b)%tile(m,n)%dx_s
            read(10,*) blk(b)%tile(m,n)%sx_s
            read(10,*) blk(b)%tile(m,n)%j_st_ism, blk(b)%tile(m,n)%j_end_ism
          end if

          read(10,*) blk(b)%tile(m,n)%j_str_smooth
          if ( blk(b)%tile(m,n)%j_str_smooth == 'yes' ) then
            read(10,*) blk(b)%tile(m,n)%dy_s
            read(10,*) blk(b)%tile(m,n)%sy_s
            read(10,*) blk(b)%tile(m,n)%s_f
          end if

          read(10,*) blk(b)%tile(m,n)%smooth_low
          read(10,*) blk(b)%tile(m,n)%smooth_upp

          read(10,*) blk(b)%tile(m,n)%connect
          if ( blk(b)%tile(m,n)%connect == 'yes' ) then
            read(10,*) blk(b)%tile(m,n)%ind_con(:)
          end if

          ! lower boundary
          read(10,*) blk(b)%tile(m,n)%np_low
          allocate( blk(b)%tile(m,n)%xp_low(blk(b)%tile(m,n)%np_low) )
          allocate( blk(b)%tile(m,n)%yp_low(blk(b)%tile(m,n)%np_low) )
          do j = 1, blk(b)%tile(m,n)%np_low
            read(10,*) blk(b)%tile(m,n)%xp_low(j), blk(b)%tile(m,n)%yp_low(j)
          end do

          ! upper boundary
          read(10,*) blk(b)%tile(m,n)%np_upp
          allocate( blk(b)%tile(m,n)%xp_upp(blk(b)%tile(m,n)%np_upp) )
          allocate( blk(b)%tile(m,n)%yp_upp(blk(b)%tile(m,n)%np_upp) )
          do j = 1, blk(b)%tile(m,n)%np_upp
            read(10,*) blk(b)%tile(m,n)%xp_upp(j), blk(b)%tile(m,n)%yp_upp(j)
          end do
        
          write(*,*) '    -> done'
        end do
      end do

    end do

  close(10)
  ! - - - - - - - - - - - - - - - - - - - - - - - - -

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: x-constant lines
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, nblocks
    do n = 1, blk(b)%n
      do m = 1, blk(b)%m

        if ( blk(b)%tile(m,n)%i_str_smooth /= 'yes' ) then

          do j = 1, blk(b)%tile(m,n)%nj+1

            blk(b)%tile(m,n)%x(:,j) = &
            stretch ( start = blk(b)%tile(m,n)%xp_low(1),   &
                      end   = blk(b)%tile(m,n)%xp_low(      &
                              blk(b)%tile(m,n)%np_low),     &
                      n     = blk(b)%tile(m,n)%ni+1,        &
                      dir   = blk(b)%tile(m,n)%sx,          &
                      delta = blk(b)%tile(m,n)%dx ) 
          end do

        elseif ( blk(b)%tile(m,n)%i_str_smooth == 'yes' ) then
            
          allocate( x_line_1( blk(b)%tile(m,n)%ni+1 ) )
          allocate( x_line_2( blk(b)%tile(m,n)%ni+1 ) )
          
          do j = 1, blk(b)%tile(m,n)%nj+1

            x_line_1 = stretch ( start = blk(b)%tile(m,n)%xp_low(1),  &
                                 end   = blk(b)%tile(m,n)%xp_low(     &
                                         blk(b)%tile(m,n)%np_low),    &
                                 n     = blk(b)%tile(m,n)%ni+1,       &
                                 dir   = blk(b)%tile(m,n)%sx,         &
                                 delta = blk(b)%tile(m,n)%dx ) 

            x_line_2 = stretch ( start = blk(b)%tile(m,n)%xp_low(1),  &
                                 end   = blk(b)%tile(m,n)%xp_low(     &
                                         blk(b)%tile(m,n)%np_low),    &
                                 n     = blk(b)%tile(m,n)%ni+1,       &
                                 dir   = blk(b)%tile(m,n)%sx_s,       &
                                 delta = blk(b)%tile(m,n)%dx_s ) 
            
            weight_1 = float ( ( j - blk(b)%tile(m,n)%j_st_ism ) ) / &
                       float ( blk(b)%tile(m,n)%j_end_ism - blk(b)%tile(m,n)%j_st_ism ) 
            weight_1 = 2d0*weight_1**3 - 3d0*weight_1**2 + 1d0
            if ( j < blk(b)%tile(m,n)%j_st_ism  ) weight_1 = 1d0
            if ( j > blk(b)%tile(m,n)%j_end_ism ) weight_1 = 0d0
            weight_2 = 1d0 - weight_1
            blk(b)%tile(m,n)%x(:,j) = x_line_1 * weight_1 + x_line_2 * weight_2

          end do

          deallocate( x_line_1, x_line_2 )

        end if
      
      end do
    end do
  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: y discretization
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, nblocks
    do n = 1, blk(b)%n
      do m = 1, blk(b)%m

        allocate( x1( blk(b)%tile(m,n)%ni+1 ) )
        allocate( y1( blk(b)%tile(m,n)%ni+1 ) )
        allocate( y2( blk(b)%tile(m,n)%ni+1 ) )

        x1 = blk(b)%tile(m,n)%x(:,1)
        call interp1( blk(b)%tile(m,n)%xp_low, blk(b)%tile(m,n)%yp_low, x1, y1 )
        x1 = blk(b)%tile(m,n)%x(:,blk(b)%tile(m,n)%nj+1)
        call interp1( blk(b)%tile(m,n)%xp_upp, blk(b)%tile(m,n)%yp_upp, x1, y2 )

        if ( blk(b)%tile(m,n)%smooth_low == 'yes' ) then
          point_1%x = x1(1)
          point_2%x = x1(size(x1))
          point_1%y = blk(b)%tile(m,n)%yp_low(1)
          point_2%y = blk(b)%tile(m,n)%yp_low(size(blk(b)%tile(m,n)%yp_low))
          call cubic( point_1%x, point_1%y, point_2%x, point_2%y, x1, y1 )
        end if

        if ( blk(b)%tile(m,n)%smooth_upp == 'yes' ) then
          point_1%x = x1(1)
          point_2%x = x1(size(x1))
          point_1%y = blk(b)%tile(m,n)%yp_upp(1)
          point_2%y = blk(b)%tile(m,n)%yp_upp(size(blk(b)%tile(m,n)%yp_upp))
          call cubic( point_1%x, point_1%y, point_2%x, point_2%y, x1, y2 )
        end if

        if ( blk(b)%tile(m,n)%j_str_smooth /= 'yes' ) then
    
          do i = 1, blk(b)%tile(m,n)%ni + 1

            blk(b)%tile(m,n)%y(i,:) =                 &
            stretch ( start = y1(i),                  &
                      end   = y2(i),                  &
                      n     = blk(b)%tile(m,n)%nj+1,  &
                      dir   = blk(b)%tile(m,n)%sy,    &
                      delta = blk(b)%tile(m,n)%dy )
          end do
    
        elseif ( blk(b)%tile(m,n)%j_str_smooth=='yes' .and. blk(b)%tile(m,n)%connect/='yes' ) then

          allocate( y_line_1( blk(b)%tile(m,n)%nj+1 ) )
          allocate( y_line_2( blk(b)%tile(m,n)%nj+1 ) )

          do i = 1, blk(b)%tile(m,n)%ni + 1
            y_line_1 = stretch ( start = y1(i),                  &
                                 end   = y2(i),                  &
                                 n     = blk(b)%tile(m,n)%nj+1,  &
                                 dir   = blk(b)%tile(m,n)%sy,    &
                                 delta = blk(b)%tile(m,n)%dy ) 

            y_line_2 = stretch ( start = y1(i),                  &
                                 end   = y2(i),                  &
                                 n     = blk(b)%tile(m,n)%nj+1,  &
                                 dir   = blk(b)%tile(m,n)%sy_s,  &
                                 delta = blk(b)%tile(m,n)%dy_s ) 

            weight_1 = ( x1(i) - x1(1) ) / ( x1(size(x1)) - x1(1) )
            weight_1 = 2d0*weight_1**3 - 3d0*weight_1**2 + 1d0
            weight_2 = 1d0 - weight_1
            blk(b)%tile(m,n)%y(i,:) = y_line_1 * weight_1 + y_line_2 * weight_2

          end do

          deallocate( y_line_1, y_line_2 )
    
        elseif ( blk(b)%tile(m,n)%j_str_smooth=='yes' .and. blk(b)%tile(m,n)%connect=='yes' ) then

          allocate( y_line_1( blk(b)%tile(m,n)%nj+1 ) )
          allocate( y_line_2( blk(b)%tile(m,n)%nj+1 ) )
          allocate( y_line_3( blk(b)%tile(m,n)%nj+1 ) )

          do i = 1, blk(b)%tile(m,n)%ni + 1
            ! Connect y_line_1 with other tile
            bs = blk(b)%tile(m,n)%ind_con(1)
            ms = blk(b)%tile(m,n)%ind_con(2)
            ns = blk(b)%tile(m,n)%ind_con(3)
            conn_nj = blk(bs)%tile(ms,ns)%nj
            conn_ni = blk(bs)%tile(ms,ns)%ni
            loca_nj = blk(b)%tile(m,n)%nj

            ! Stretched line representative of the connected block
            y_line_1(1:conn_nj+1) = blk(bs)%tile(ms,ns)%y(conn_ni+1,:)
            y_line_3(1:conn_nj+1) = linspace(y_line_1(1), y_line_1(conn_nj+1), conn_nj+1)
            weight_1 = ( x1(i) - x1(1) ) / ( x1(size(x1)) - x1(1) )
            weight_1 = 2d0*weight_1**3 - 3d0*weight_1**2 + 1d0
            weight_1 = max( weight_1, 0d0 )
            weight_2 = 1d0 - weight_1
            !y_line_1(1:conn_nj+1) = y_line_1(1:conn_nj+1)*weight_1 + &
            !                        y_line_3(1:conn_nj+1)*weight_2
            
            ! Damp y_line_1 to take into account y_upp of local wall
            y_line_1(1:conn_nj+1) = ( y_line_1(1:conn_nj+1) - y_line_1(1) ) &
                                  * ( y2(i) - y1(i) )/( y2(1) - y1(1) )     &
                                  + y_line_1(1)


            y_line_1(conn_nj+1:loca_nj+1) =         &
            stretch ( start = y_line_1(conn_nj+1),  &
                      end   = y2(i),                &
                      n     = loca_nj-conn_nj+1,    &
                      dir   = blk(b)%tile(m,n)%sy,  &
                      delta = blk(b)%tile(m,n)%dy ) 

            y_line_2 = stretch ( start = y1(i),                  &
                                 end   = y2(i),                  &
                                 n     = blk(b)%tile(m,n)%nj+1,  &
                                 dir   = blk(b)%tile(m,n)%sy_s,  &
                                 delta = blk(b)%tile(m,n)%dy_s ) 

            blk(b)%tile(m,n)%y(i,:) = y_line_1 * weight_1 + y_line_2 * weight_2

          end do

          deallocate( y_line_1, y_line_2, y_line_3 )

        end if

    
        deallocate( x1, y1, y2 )

      end do
    end do
  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: adjust j-lines
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, nblocks
    do n = 1, blk(b)%n
      do m = 1, blk(b)%m

        trapez_1 = blk(b)%tile(m,n)%xp_low(1) /= blk(b)%tile(m,n)%xp_upp(1)
        trapez_2 = blk(b)%tile(m,n)%xp_low(blk(b)%tile(m,n)%np_low) /= &
                   blk(b)%tile(m,n)%xp_upp(blk(b)%tile(m,n)%np_upp)

        if ( .not.trapez_1 .and. .not. trapez_2 ) cycle

        lenght_low = blk(b)%tile(m,n)%xp_low(blk(b)%tile(m,n)%np_low) &
                   - blk(b)%tile(m,n)%xp_low(1)
        lenght_upp = blk(b)%tile(m,n)%xp_upp(blk(b)%tile(m,n)%np_upp) &
                   - blk(b)%tile(m,n)%xp_upp(1)

        do j = 1, blk(b)%tile(m,n)%nj + 1

          weight_1 = blk(b)%tile(m,n)%y(1,blk(b)%tile(m,n)%nj+1) - blk(b)%tile(m,n)%y(1,1)
          weight_1 = ( blk(b)%tile(m,n)%y(1,j) - blk(b)%tile(m,n)%y(1,1) ) / weight_1
          weight_2 = 1d0 - weight_1
          lenght_loc = weight_2 * lenght_low + weight_1 * lenght_upp     
          blk(b)%tile(m,n)%x(:,j) = ( blk(b)%tile(m,n)%x(:,j) - blk(b)%tile(m,n)%x(1,j) ) &
                                  / ( blk(b)%tile(m,n)%x(blk(b)%tile(m,n)%ni+1,j) &
                                  -   blk(b)%tile(m,n)%x(1,j) ) &
                                  * lenght_loc + blk(b)%tile(m,n)%x(1,j)
        
        end do

      end do
    end do
  end do


  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Block dimensions
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  write(*,*), ' Assembling blocks'
  do b = 1, nblocks
    blk(b)%ni = 0
    blk(b)%nj = 0
    n = 1
    do m = 1, blk(b)%m
      blk(b)%ni = blk(b)%ni + blk(b)%tile(m,n)%ni
    end do
    m = 1
    do n = 1, blk(b)%n
      blk(b)%nj = blk(b)%nj + blk(b)%tile(m,n)%nj
    end do
    write(*,*) '  Block ', b, blk(b)%ni,'x',blk(b)%nj
    allocate( blk(b)%x(blk(b)%ni+1,blk(b)%nj+1) )
    allocate( blk(b)%y(blk(b)%ni+1,blk(b)%nj+1) )
  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Block assembling
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, nblocks

    cursor_j = 1
    do n = 1, blk(b)%n

      cursor_i = 1
      do m = 1, blk(b)%m

        blk(b)%x( cursor_i : cursor_i + blk(b)%tile(m,n)%ni,   &
                  cursor_j : cursor_j + blk(b)%tile(m,n)%nj) = &
        blk(b)%tile(m,n)%x

        blk(b)%y( cursor_i : cursor_i + blk(b)%tile(m,n)%ni,   &
                  cursor_j : cursor_j + blk(b)%tile(m,n)%nj) = &
        blk(b)%tile(m,n)%y

        cursor_i = cursor_i + blk(b)%tile(m,n)%ni

      end do

      cursor_j = cursor_j + blk(b)%tile(m-1,n)%nj

    end do

  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Writing the mesh
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  open(10,file='mesh.dat',status='unknown')
  write(10,*) 'TITLE     = "Mesh"'
  write(10,*) 'VARIABLES = "x"'
  write(10,*) '"y"'

  do b = 1, nblocks

    write(10,*)'ZONE T="Block ', b,'"'
    
    write(10,mesh_dim) blk(b)%ni+1, blk(b)%nj+1
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'DT=(SINGLE SINGLE)'

    do j = 1, blk(b)%nj + 1
      do i = 1, blk(b)%ni + 1
        write (10,real_form) blk(b)%x(i,j)
      end do
    end do

    do j = 1, blk(b)%nj + 1
      do i = 1, blk(b)%ni +1
        write (10,real_form) blk(b)%y(i,j)
      end do
    end do

  end do
  close(10)

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Total number of cells
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  ncell = 0
  do b = 1, nblocks
    ncell = ncell + blk(b)%ni*blk(b)%nj
  end do
  write(*,*) ' Total number of cells', ncell

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Writing the coarse meshes
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  write(*,*) ' 2nd level mesh with: ', ncell/4, 'cells'
  open(10,file='mesh2.dat',status='unknown')
  write(10,*) 'TITLE     = "Mesh"'
  write(10,*) 'VARIABLES = "x"'
  write(10,*) '"y"'

  do b = 1, nblocks

    write(10,*)'ZONE T="Block ', b,'"'
    
    write(10,mesh_dim) blk(b)%ni/2+1, blk(b)%nj/2+1
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'DT=(SINGLE SINGLE)'

    do j = 1, blk(b)%nj + 1, 2
      do i = 1, blk(b)%ni + 1, 2
        write (10,real_form) blk(b)%x(i,j)
      end do
    end do

    do j = 1, blk(b)%nj + 1, 2
      do i = 1, blk(b)%ni + 1, 2
        write (10,real_form) blk(b)%y(i,j)
      end do
    end do

  end do
  close(10)

  write(*,*) ' 3nd level mesh with: ', ncell/16, 'cells'
  open(10,file='mesh3.dat',status='unknown')
  write(10,*) 'TITLE     = "Mesh"'
  write(10,*) 'VARIABLES = "x"'
  write(10,*) '"y"'

  do b = 1, nblocks

    write(10,*)'ZONE T="Block ', b,'"'
    
    write(10,mesh_dim) blk(b)%ni/4+1, blk(b)%nj/4+1
    write(10,*) 'DATAPACKING=BLOCK'
    write(10,*) 'DT=(SINGLE SINGLE)'

    do j = 1, blk(b)%nj + 1, 4
      do i = 1, blk(b)%ni + 1, 4
        write (10,real_form) blk(b)%x(i,j)
      end do
    end do

    do j = 1, blk(b)%nj + 1, 4
      do i = 1, blk(b)%ni + 1, 4
        write (10,real_form) blk(b)%y(i,j)
      end do
    end do

  end do
  close(10)

end program patchwork