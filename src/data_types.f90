module Data_Types
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

  type :: Tile_Type
    integer :: ni, nj, np_low, np_upp, ind_con(3)
    integer :: j_st_ism, j_end_ism
    character(2) :: sx, sy
    character(3) :: i_str_smooth, j_str_smooth, smooth_low, smooth_upp, connect
    character(2) :: sy_s, sx_s
    real(R8) :: dx, dy, dx_s, dy_s, s_f
    real(R8), allocatable :: xp_low(:), yp_low(:)
    real(R8), allocatable :: xp_upp(:), yp_upp(:)
    real(R8), allocatable :: x(:,:), y(:,:)
  end type Tile_Type

  type :: Point_Type
    real(R8) :: x, y
  end type Point_Type

  type :: Block_Type
    integer :: m, n
    integer :: ni, nj
    real(R8), allocatable :: x(:,:), y(:,:)
    type(Tile_Type), allocatable :: tile(:,:)
  end type

  type :: Grid_Type
    integer :: nblocks
    type(Block_Type), allocatable :: blk(:)
  end type

 end module Data_Types