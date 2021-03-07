program wave_function

    implicit none

    integer, parameter :: dp = selected_real_kind(14,200)
    real(dp), parameter :: hbar = 1_dp, m_e = 1_dp
    real(dp) :: E, x_beg, x_end, delta_h, dirac_delta_pot
    integer(dp) :: grid, i, l
    !real :: start, finish
    
    
    real(dp), allocatable :: x(:) , u_wave(:), f(:), y(:)

    grid = 10000
    E = 1 ! Energy of the system in Atomic Units
    l = 10 ! l value for the effective potential
    x_beg = 1d-40 ! Begining of spataial coordinates
    x_end = 50_dp ! End of spatial coordiantes
    delta_h = (x_end-x_beg)/grid

    allocate(x(0:grid), u_wave(0:grid), f(0:grid), y(0:grid))


    ! Loop - 1 : calculating the function potential for the 
    ! schrodinger equation.

    x(0) = x_beg
    
    do i = 0, grid
        f(i) = -(2*m_e/hbar**2)*(E- dirac_delta_pot(x(i))-((hbar/(2*m_e))*(l*(l+1)/x(i)**2)))
        y(i) = 1 - (delta_h**2)/(12)*f(i)
        x(i+1) = x(i) + delta_h
        !print '(3e16.8, 3e16.8)',x(i), y(i)
    end do

    ! Loop - 2 Running Schroedinger's equation using
    ! Numerov's method.

    u_wave(0) = x(0)**(l+1)
    u_wave(1) = x(1)**(l+1)

    do i = 1, grid
        u_wave(i+1) = (u_wave(i)*(12.0_dp-(10.0_dp*y(i)))-y(i-1)*u_wave(i-1))/y(i+1)
    end do

    !u_wave(grid/2) = 1d-21
    !u_wave(grid/2+1) = delta_h

    !do i = grid/2+1, grid
        !u_wave(i+1) = (u_wave(i)*(12.0_dp-(10.0_dp*y(i)))-y(i-1)*u_wave(i-1))/y(i+1)
    !end do

    u_wave=u_wave/sqrt(dot_product(u_wave,u_wave)*delta_h)

    
    do  i = 0, grid
        print '(3e16.8, 3e16.8)',x(i), u_wave(i)
    end do

end program wave_function

real*8 function dirac_delta_pot(x)

    implicit none
    integer, parameter :: dp = selected_real_kind(14,200)
    real(dp), intent(in) :: x
    real(dp) :: a, pi, x_end
    x_end = 50_dp
    a = 1d-40
    pi = 22_dp/7_dp
    dirac_delta_pot = (1_dp/(a*sqrt(pi)))*exp(-((x-(x_end/2))/a)**2)

end function dirac_delta_pot
