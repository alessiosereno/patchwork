module mesh_writer_mod
  use Data_Types
  use IO_mod, only: dl
  implicit none
  private
  public :: write_mesh, write_tiles

  character(8),  parameter :: real_form = '(e20.10)'
  character(49), parameter :: mesh_dim  = "(' I= ',i3,', J= ',i3,', K= 1, ZONETYPE=Ordered')"

contains

  ! Writes a Tecplot-ordered mesh with optional subsampling stride.
  ! stride=1 -> full mesh ; stride=2 or 4 -> coarser levels.
  subroutine write_mesh(grid, filename, stride)
    type(Grid_Type),   intent(in) :: grid
    character(len=*),  intent(in) :: filename
    integer,           intent(in) :: stride

    integer :: b, i, j, ni_out, nj_out

    open(10, file=filename, status='unknown')
    write(10,*) 'TITLE     = "Mesh"'
    write(10,*) 'VARIABLES = "x"'
    write(10,*) '"y"'

    do b = 1, grid%nblocks
      ni_out = grid%blk(b)%ni/stride + 1
      nj_out = grid%blk(b)%nj/stride + 1

      write(10,*) 'ZONE T="Block ', b, '"'
      write(10, mesh_dim) ni_out, nj_out
      write(10,*) 'DATAPACKING=BLOCK'
      write(10,*) 'DT=(SINGLE SINGLE)'

      do j = 1, grid%blk(b)%nj + 1, stride
        do i = 1, grid%blk(b)%ni + 1, stride
          write(10, real_form) grid%blk(b)%x(i,j)
        end do
      end do
      do j = 1, grid%blk(b)%nj + 1, stride
        do i = 1, grid%blk(b)%ni + 1, stride
          write(10, real_form) grid%blk(b)%y(i,j)
        end do
      end do
    end do
    close(10)
  end subroutine write_mesh


  subroutine write_tiles(grid, filename)
    type(Grid_Type),   intent(in) :: grid
    character(len=*),  intent(in) :: filename

    integer :: b, m, n, i, j

    open(10, file=filename, status='unknown')
    write(10,*) 'TITLE     = "Tiles"'
    write(10,*) 'VARIABLES = "x"'
    write(10,*) '"y"'

    do b = 1, grid%nblocks
      do m = 1, grid%blk(b)%m
        do n = 1, grid%blk(b)%n
          associate( t => grid%blk(b)%tile(m,n) )
          write(10,'(A16,I3,A6,I3,I3,A)') ' ZONE T= "Block:', b, ' Tile:', m, n, '"'
          write(10, mesh_dim) t%ni+1, t%nj+1
          write(10,*) 'DATAPACKING=BLOCK'
          write(10,*) 'DT=(SINGLE SINGLE)'

          do j = 1, t%nj+1
            do i = 1, t%ni+1
              write(10, real_form) t%x(i,j)
            end do
          end do
          do j = 1, t%nj+1
            do i = 1, t%ni+1
              write(10, real_form) t%y(i,j)
            end do
          end do
          end associate
        end do
      end do
    end do
    close(10)
  end subroutine write_tiles

end module mesh_writer_mod
