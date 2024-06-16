module stretching_mod
    use iso_fortran_env, only: dp => real64
    use linspace_mod
    implicit none
contains

    ! TANH stretching functions
    function stretch ( Type, Par, N, xsi, delta, start, end, dir ) result ( discr )
      ! Stretching tangente iperbolica. 
      ! Input: tipo stretching, parametro di stretching e numero di punti
      ! Out: vettore da 0 a 1 stretchato.    
      implicit none
      character, optional :: Type
      character(2), optional :: dir
      integer, intent(in) :: n
      real(dp), optional :: Par, delta, xsi, start, end
      real(dp) :: discr(N)
      ! Local
      character :: sType
      character(2) :: sDir
      integer, parameter :: try_max = 100
      integer :: try
      real(dp), parameter :: tol = 1d-3, zero = 0d0, uno = 1d0, half = 0.5d0, due = 2d0
      real(dp) :: left, right, norm, delta_n, A, xs, x(N), delta_try, error, A_1, A_2, A_try, discr_try(n)

      ! Setting
      if ( .not. present( start ) .and. .not. present( end ) ) then
        left = zero
        right = uno
      else
        left = start
        right = end
      end if
      norm = right - left

      if ( .not. present( Par ) ) then
        A = zero
      else
        A = Par
      end if

      if ( .not. present( xsi ) ) then
        xs = half
      else
        xs = xsi
      end if 

      if ( .not. present( dir ) ) then
        sDir = '??'
      else
        sDir = dir
      end if

      if ( .not. present( Type ) )  then
        sType = 'N'
      else
        sType = Type
      end if

      x = linspace ( zero, uno, n )

      if ( present( delta ) ) then

        delta_n = delta / norm

        A_1 = 0.1d0
        A_2 = 8d0
        A_try = 0.2d0

        discr_try = Tanh( A_try * x ) / Tanh( A_try )
        if ( sDir == '<>' ) &
          discr_try = ( Tanh( A_try * ( x - xs ) ) + Tanh( A_try * xs ) ) / ( Tanh( A_try * ( uno - xs ) ) + Tanh( A_try * xs ) )
        delta_try = discr_try(n) - discr_try(n-1)
        
        try = 0
        error = Abs( delta_try -  delta_n ) / delta_n

        do while ( error > tol .and. try < try_max )

          if ( delta_try < delta_n ) then
            A_2 = A_try
          elseif ( delta_try > delta_n ) then
            A_1 = A_try
          end if

          A_try = 0.5d0 * ( A_1 + A_2 )
          discr_try = Tanh( A_try * x ) / Tanh( A_try )
          if ( sDir == '<>' ) &
            discr_try = ( Tanh( A_try * ( x - xs ) ) + Tanh( A_try * xs ) ) / ( Tanh( A_try * ( uno - xs ) ) + Tanh( A_try * xs ) )
          delta_try = discr_try(n) - discr_try(n-1)
          error = Abs( delta_try -  delta_n ) / delta_n
          try = try + 1
        end do

        A = A_try

        if ( sDir == '<-' ) A = - A

      end if


      if ( sType == 'A' .or. sDir == '->' .or. sDir == '<-' ) then

        if ( A > zero .or. dir == '->') then
          discr = Tanh( A * x ) / Tanh( A )

        elseif ( A < zero .or. sDir == '<-') then
          discr = uno + Tanh( Abs(A) * ( x - uno ) ) / Tanh( Abs(A) )

        end if

      elseif ( sType == 'S' .or. sDir == '<>' .or. sDir == '><' .or. sDir == '==' ) then

        if ( A > zero .or. sDir == '<>' ) then
          discr = ( Tanh( A * ( x - xs ) ) + Tanh( A * xs ) ) / ( Tanh( A * ( uno - xs ) ) + Tanh( A * xs ) )

        elseif ( A < zero .or. sDir == '><' ) then 
          discr = due * x - ( Tanh( Abs(A) * ( x - xs ) ) + Tanh ( Abs(A) * xs) ) / ( Tanh( Abs(A) * ( uno - xs ) ) + Tanh( Abs(A) * xs ) )
        
        elseif ( sDir == '==' ) then
          discr = x

        end if

      end if

      discr = discr * norm + left

    end function stretch

end module stretching_mod
