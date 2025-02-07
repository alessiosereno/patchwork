module Tec3D_mod
    use finer, only : file_ini
    use iso_fortran_env, only: I4 => int32, R8 => real64
    use Data_Types
    use linspace_mod
    implicit none
    contains
  
    subroutine Tec3D (grid)
      implicit none
      type(Grid_Type), intent(in) :: grid
      integer :: Nb, Nx, Ny, Nz, i, j, k, b, Nz_fix, error
      integer, allocatable :: Nx_arr(:), Ny_arr(:), Nz_arr(:)
      real(R8), allocatable :: X(:,:,:,:), Y(:,:,:,:), Z(:,:,:,:)
      real(R8), allocatable :: a(:), zdir(:)
      real(R8) :: axis, zmax, c, p, theta_max, width_max
      character(len=256) :: inputfile, outputfile, mesh_type, stretch
      logical :: threedim, circular_sector
      type(file_ini) :: fini
  
      inputfile = 'mesh.dat'
      outputfile = 'mesh.tec'
      circular_sector = .false.
      threedim = .false.
      axis = 0_R8

      call fini%load(filename='input.ini')
      call fini%get(section_name='PATCHWORK', option_name='inputfile', val=inputfile, error=error)
      call fini%get(section_name='PATCHWORK', option_name='outputfile', val=outputfile, error=error)
      call fini%get(section_name='PATCHWORK', option_name='type', val=mesh_type, error=error)
      if (error/=0) mesh_type = '3D'
      call fini%get(section_name='PATCHWORK', option_name='axis', val=axis, error=error)

      select case(mesh_type)
      case('2D')
        Nz_fix = 2
        call fini%get(section_name='PATCHWORK', option_name='width', val=width_max, error=error)
        if (error/=0) then
          width_max = 1e-3
          zmax = width_max
          stretch = 'none'
        end if
      case('2Daxi')
        Nz_fix = 2
        call fini%get(section_name='PATCHWORK', option_name='theta', val=theta_max, error=error)
        if (error/=0) then
          theta_max = 3.141592653589793_R8 / 180_R8
          zmax = theta_max
          circular_sector = .true.
        end if
        stretch = 'none'
      case('3D')
        call fini%get(section_name='PATCHWORK', option_name='Nz', val=Nz_fix, error=error)
        if (error==0) then
          Nz_fix = Nz_fix + 1
          call fini%get(section_name='PATCHWORK', option_name='theta', val=theta_max, error=error)
          if (error==0) then
            theta_max = theta_max * 3.141592653589793_R8 / 180_R8
            zmax = theta_max
            circular_sector = .true.
          else
            call fini%get(section_name='PATCHWORK', option_name='width', val=width_max, error=error)
            if (error/=0) error stop ( '  ERROR : theta or width parameter necessary!' )
            zmax = width_max
          end if
          call fini%get(section_name='PATCHWORK', option_name='stretch', val=stretch, error=error)
          if (error==0) then
            call fini%get(section_name='PATCHWORK', option_name='c', val=c, error=error)
            call fini%get(section_name='PATCHWORK', option_name='p', val=p, error=error)
          else
            stretch = 'none'
          end if
        else
          threedim = .true.
          write(*,*) '  Look for a 3D mesh file'
        end if
      end select

      if (.not.threedim) then
        a = linspace(0d0,real(Nz_fix,kind=8),Nz_fix)/Nz_fix
        if (trim(stretch) == 'none') then
          ! No stretch
          zdir = zmax*a-0.5_R8*zmax
          print*, 'no stretch ok'
        elseif (trim(stretch) == 'left') then
          ! Left
          zdir = zmax*tanh(c*a)/tanh(c)-0.5*zmax
        elseif (trim(stretch) == 'outer') then
          ! Outer
          zdir = zmax*(tanh(c*(a-p))+tanh(c*p))/(tanh(c*(1-p))+tanh(c*p))-0.5*zmax
        end if
      end if

      ! 3D mesh
      Nb = grid%nblocks
     
      allocate(Nx_arr(Nb), Ny_arr(Nb), Nz_arr(Nb))
     
      ! Read block dimensions
      do b = 1, Nb
        Nx_arr(b) = grid%blk(b)%ni+1
        Ny_arr(b) = grid%blk(b)%nj+1
        Nz_arr(b) = Nz_fix
      end do
     
      ! Find max dimensions to allocate arrays
      Nx = maxval(Nx_arr)
      Ny = maxval(Ny_arr)
      Nz = maxval(Nz_arr)
     
      allocate(X(Nx,Ny,Nz,Nb), Y(Nx,Ny,Nz,Nb), Z(Nx,Ny,Nz,Nb))
     
      ! Read mesh data
      do b = 1, Nb
        do j = 1, Ny_arr(b)
          do i = 1, Nx_arr(b)
            X(i,j,1,b) = grid%blk(b)%x(i,j)
          end do
        end do
        do j = 1, Ny_arr(b)
          do i = 1, Nx_arr(b)
            Y(i,j,1,b) = grid%blk(b)%y(i,j)
            if (Y(i,j,1,b) <= axis) Y(i,j,1,b) = axis
          end do
        end do
      end do
     
      ! Handle circular sector transformation
      if (circular_sector) then
        do b = 1, Nb
          do j = 1, Ny_arr(b)
            do i = 1, Nx_arr(b)
              Z(i,j,1,b) = Y(i,j,1,b) !* sin(theta_max)
              Y(i,j,1,b) = Y(i,j,1,b) !* cos(theta_max)
            end do
          end do
        end do
      end if
     
      ! ! Stretching in z-direction
      ! do b = 1, Nb
      !   do k = 1, Nz_arr(b)
      !     Z(:,:,k,b) = (real(k-1)/real(Nz_arr(b)-1)) * zmax - 0.5_R8 * zmax
      !   end do
      ! end do
     
      ! Open output file
      open(20, file=outputfile, status='unknown', action='write')
      write(20, '(A)') 'TITLE = "PATCHWORK mesh"'
      write(20, '(A)') 'VARIABLES = "X", "Y", "Z"'
     
      ! Write mesh data
      do b = 1, Nb
        write(20, '(A,I4,A,I4,A,I4,A)') 'ZONE I=', Nx_arr(b), ', J=', Ny_arr(b), ', K=', Nz_arr(b), ', ZONETYPE=Ordered'
        write(20, '(A)') 'DATAPACKING=BLOCK'
       
        do k = 1, Nz_arr(b)
          do j = 1, Ny_arr(b)
            do i = 1, Nx_arr(b)
              write(20, '(F12.6)') X(i,j,k,b)
            end do
          end do
        end do
       
        do k = 1, Nz_arr(b)
          do j = 1, Ny_arr(b)
            do i = 1, Nx_arr(b)
              write(20, '(F12.6)') Y(i,j,k,b)
            end do
          end do
        end do
       
        do k = 1, Nz_arr(b)
          do j = 1, Ny_arr(b)
            do i = 1, Nx_arr(b)
              write(20, '(F12.6)') Z(i,j,k,b)
            end do
          end do
        end do
      end do
      close(20)
     
      ! Cleanup
      deallocate(X, Y, Z, Nx_arr, Ny_arr, Nz_arr)
     
      print *, 'Mesh processing complete.'
    end subroutine Tec3D
  
  end module Tec3D_mod