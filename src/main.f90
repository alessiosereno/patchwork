program patchwork
  use Data_Types
  use finer, only: file_ini
  use IO_mod
  use linspace_mod
  use stretching_mod
  use interp1_mod
  use spline_mod
  use cubic_mod
  use Tec3D_mod
  use input_mod
  use mesh_writer_mod
  implicit none

  logical :: trapez_1, trapez_2
  integer :: b, i, j, m, n, bs, ncell
  integer :: js1, js2, cursor_i, cursor_j
  character(len=dl) :: OUT_Path
  logical :: write_coarse
  real(R8), allocatable, dimension(:) :: x1, y1, y2, y_line_1, y_line_2, x_line_1, x_line_2
  real(R8) :: weight_1, weight_2
  real(R8) :: lenght_loc, lenght_low, lenght_upp
  integer  :: conn_nj, conn_ni, loca_nj
  type(Grid_Type) :: grid
  type(Point_Type) :: point_1, point_2
  type(file_ini) :: fini

  call read_input(grid, fini, OUT_Path, write_coarse)

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Block dimensions
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  write(*,*) ' Allocating blocks.'
  do b = 1, grid%nblocks
    associate( blk => grid%blk(b) )
    blk%ni = 0
    blk%nj = 0
    do m = 1, blk%m
      blk%ni = blk%ni + blk%tile(m,1)%ni
    end do
    do n = 1, blk%n
      blk%nj = blk%nj + blk%tile(1,n)%nj
    end do
    write(*,*) '  Block ', b, blk%ni, 'x', blk%nj
    allocate( blk%x(blk%ni+1, blk%nj+1) )
    allocate( blk%y(blk%ni+1, blk%nj+1) )
    end associate
  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: x-constant lines
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, grid%nblocks
    do n = 1, grid%blk(b)%n
      do m = 1, grid%blk(b)%m
        associate( t => grid%blk(b)%tile(m,n) )

        if ( t%i_str_smooth /= 'yes' ) then

          do j = 1, t%nj+1
            t%x(:,j) = stretch( start = t%xp_low(1),       &
                                end   = t%xp_low(t%np_low), &
                                n     = t%ni+1,             &
                                dir   = t%sx,               &
                                delta = t%dx )
          end do

        else

          allocate( x_line_1(t%ni+1), x_line_2(t%ni+1) )

          do j = 1, t%nj+1
            x_line_1 = stretch( start = t%xp_low(1),       &
                                end   = t%xp_low(t%np_low), &
                                n     = t%ni+1,             &
                                dir   = t%sx,               &
                                delta = t%dx )
            x_line_2 = stretch( start = t%xp_low(1),       &
                                end   = t%xp_low(t%np_low), &
                                n     = t%ni+1,             &
                                dir   = t%sx_s,             &
                                delta = t%dx_s )

            weight_1 = float( j - t%j_st_ism ) / float( t%j_end_ism - t%j_st_ism )
            weight_1 = 2d0*weight_1**3 - 3d0*weight_1**2 + 1d0
            if ( j < t%j_st_ism  ) weight_1 = 1d0
            if ( j > t%j_end_ism ) weight_1 = 0d0
            weight_2 = 1d0 - weight_1
            t%x(:,j) = x_line_1 * weight_1 + x_line_2 * weight_2
          end do

          deallocate( x_line_1, x_line_2 )

        end if

        end associate
      end do
    end do
  end do

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Tile mesh building: y discretization
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  do b = 1, grid%nblocks

    grid%blk(b)%ni = 0
    grid%blk(b)%nj = 0

    do n = 1, grid%blk(b)%n
      do m = 1, grid%blk(b)%m
        associate( t => grid%blk(b)%tile(m,n) )

        allocate( x1(t%ni+1), y1(t%ni+1), y2(t%ni+1) )

        ! Interpolate lower and upper boundaries of the tile
        x1 = t%x(:,1)
        if ( t%np_low < 10 ) then
          call interp1( t%xp_low, t%yp_low, x1, y1 )
        else
          y1 = spline( t%xp_low, t%yp_low, x1 )
        end if
        x1 = t%x(:, t%nj+1)
        if ( t%np_upp < 10 ) then
          call interp1( t%xp_upp, t%yp_upp, x1, y2 )
        else
          y2 = spline( t%xp_upp, t%yp_upp, x1 )
        end if

        ! Cubic transition option on lower/upper segments
        if ( t%smooth_low == 'yes' ) then
          point_1%x = x1(1)
          point_2%x = x1(size(x1))
          point_1%y = t%yp_low(1)
          point_2%y = t%yp_low(size(t%yp_low))
          call cubic( point_1%x, point_1%y, point_2%x, point_2%y, x1, y1 )
        end if
        if ( t%smooth_upp == 'yes' ) then
          point_1%x = x1(1)
          point_2%x = x1(size(x1))
          point_1%y = t%yp_upp(1)
          point_2%y = t%yp_upp(size(t%yp_upp))
          call cubic( point_1%x, point_1%y, point_2%x, point_2%y, x1, y2 )
        end if

        if ( t%j_str_smooth /= 'yes' .and. t%connect /= 'yes' ) then

          do i = 1, t%ni+1
            t%y(i,:) = stretch( start = y1(i),   &
                                end   = y2(i),   &
                                n     = t%nj+1,  &
                                dir   = t%sy,    &
                                delta = t%dy )
          end do

        elseif ( t%j_str_smooth == 'yes' .and. t%connect /= 'yes' ) then

          allocate( y_line_1(t%nj+1), y_line_2(t%nj+1) )

          do i = 1, t%ni+1
            y_line_1 = stretch( start = y1(i),   &
                                end   = y2(i),   &
                                n     = t%nj+1,  &
                                dir   = t%sy,    &
                                delta = t%dy )
            y_line_2 = stretch( start = y1(i),   &
                                end   = y2(i),   &
                                n     = t%nj+1,  &
                                dir   = t%sy_s,  &
                                delta = t%dy_s )

            weight_1 = ( x1(i) - x1(1) ) / ( x1(size(x1)) - x1(1) )
            if ( t%def_jstr_sm_range ) then
              weight_1 = ( x1(i) - t%x_jstr_sm_range(1) ) &
                       / ( t%x_jstr_sm_range(2) - t%x_jstr_sm_range(1) )
            end if
            if ( weight_1 < 0.1_R8 ) then
              t%y(i,:) = y_line_1
            elseif ( weight_1 > 0.9_R8 ) then
              t%y(i,:) = y_line_2
            else
              weight_1 = (2._R8*weight_1**3 - 3._R8*weight_1**2 + 0.54_R8*weight_1)/0.512_R8 + 0.949219_R8
              weight_2 = 1._R8 - weight_1
              t%y(i,:) = y_line_1 * weight_1 + y_line_2 * weight_2
            end if
          end do

          deallocate( y_line_1, y_line_2 )

        elseif ( t%connect == 'yes' ) then

          allocate( y_line_1(t%nj+1), y_line_2(t%nj+1) )

          bs      = t%ind_con(1)
          js1     = t%ind_con(2)
          js2     = t%ind_con(3)
          conn_nj = grid%blk(bs)%nj
          conn_ni = grid%blk(bs)%ni
          loca_nj = t%nj

          do i = 1, t%ni+1

            ! Discretized line inherited from connected block
            y_line_1(1:loca_nj+1) = grid%blk(bs)%y(conn_ni+1, js1:js2)
            weight_1 = ( x1(i) - x1(1) ) / ( x1(size(x1)) - x1(1) )
            weight_1 = 2._R8*weight_1**3 - 3._R8*weight_1**2 + 1._R8
            weight_1 = max( weight_1, 0._R8 )
            weight_2 = 1._R8 - weight_1

            ! Rescale inherited line onto local wall span
            y_line_1(1:loca_nj+1) = ( y_line_1(1:loca_nj+1) - y_line_1(1) )   &
                                  * ( y2(i) - y1(i) ) / ( y2(1) - y1(1) )    &
                                  + y_line_1(1)

            y_line_2 = stretch( start = y1(i),   &
                                end   = y2(i),   &
                                n     = t%nj+1,  &
                                dir   = t%sy,    &
                                delta = t%dy )

            t%y(i,:) = y_line_1 * weight_1 + y_line_2 * weight_2
          end do

          deallocate( y_line_1, y_line_2 )

        end if

        deallocate( x1, y1, y2 )

        ! Trapezoidal-shaped tile x-adjustment
        trapez_1 = t%xp_low(1)         /= t%xp_upp(1)
        trapez_2 = t%xp_low(t%np_low)  /= t%xp_upp(t%np_upp)

        if ( trapez_1 .or. trapez_2 ) then

          lenght_low = t%xp_low(t%np_low) - t%xp_low(1)
          lenght_upp = t%xp_upp(t%np_upp) - t%xp_upp(1)

          do j = 1, t%nj+1

            weight_1 = t%y(1, t%nj+1) - t%y(1,1)
            weight_1 = ( t%y(1,j) - t%y(1,1) ) / weight_1
            weight_2 = 1d0 - weight_1
            lenght_loc = weight_2 * lenght_low + weight_1 * lenght_upp

            if ( trapez_1 .and. .not. trapez_2 ) then
              t%x(:,j) = ( t%x(:,j) - t%x(1,j) )                &
                       / ( t%x(t%ni+1,j) - t%x(1,j) )           &
                       * lenght_loc + t%x(1,j)                  &
                       + (lenght_low - lenght_upp) * weight_1

            elseif ( trapez_2 .and. .not. trapez_1 ) then
              t%x(:,j) = ( t%x(:,j) - t%x(1,j) )                &
                       / ( t%x(t%ni+1,j) - t%x(1,j) )           &
                       * lenght_loc + t%x(1,j)

            elseif ( trapez_1 .and. trapez_2 ) then
              t%x(:,j) = ( t%x(:,j) - t%x(1,j) )                &
                       / ( t%x(t%ni+1,j) - t%x(1,j) )           &
                       * lenght_loc + t%x(1,j)                  &
                       + (lenght_low - lenght_upp) * weight_1 /2d0

            end if

          end do

        end if

        if ( n == 1 ) then
          grid%blk(b)%ni = grid%blk(b)%ni + t%ni
          write(*,*) ' Updating i-dimension of block', b, '->', grid%blk(b)%ni
        end if

        end associate
      end do ! m loop

      grid%blk(b)%nj = grid%blk(b)%nj + grid%blk(b)%tile(1,n)%nj
      write(*,*) ' Updating j-dimension of block', b, '->', grid%blk(b)%nj

    end do ! n loop

    write(*,*) ' Assembling block', b
    cursor_j = 1
    do n = 1, grid%blk(b)%n
      cursor_i = 1
      do m = 1, grid%blk(b)%m
        associate( t => grid%blk(b)%tile(m,n) )
        grid%blk(b)%x( cursor_i:cursor_i+t%ni, cursor_j:cursor_j+t%nj ) = t%x
        grid%blk(b)%y( cursor_i:cursor_i+t%ni, cursor_j:cursor_j+t%nj ) = t%y
        cursor_i = cursor_i + t%ni
        end associate
      end do
      cursor_j = cursor_j + grid%blk(b)%tile(m-1,n)%nj
    end do
    write(*,*) '    -> done'

  end do ! block loop

  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  !    Writing the meshes
  ! - - - - - - - - - - - - - - - - - - - - - - - - -
  call write_mesh(grid, trim(OUT_Path)//'mesh.dat', stride=1)

  ncell = 0
  do b = 1, grid%nblocks
    ncell = ncell + grid%blk(b)%ni * grid%blk(b)%nj
  end do
  write(*,*) ' Total number of cells', ncell

  if (write_coarse) then
    write(*,*) ' 2nd level mesh with: ', ncell/4, 'cells'
    call write_mesh(grid, trim(OUT_Path)//'mesh2.dat', stride=2)
    write(*,*) ' 3nd level mesh with: ', ncell/16, 'cells'
    call write_mesh(grid, trim(OUT_Path)//'mesh3.dat', stride=4)
  end if

  call write_tiles(grid, trim(OUT_Path)//'tiles.dat')

  !call Tec3D (grid) TODO fix

end program patchwork
