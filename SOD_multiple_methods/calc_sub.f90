module calc_sub
  use variables
  use ieee_arithmetic

  implicit none
contains

  ! Subroutine to initialize the grid
  subroutine grid()
    integer :: i
    dx = (x_rbound - x_lbound) / nx
    ! Initialize centers
    do i = 0, nx-1
        x_center(i+1) = x_lbound + (real(i, dp) + 0.5d0) * dx
    end do
    ! Initialize edges
    do i = 0, nx
        x_edges(i+1) = x_lbound + real(i, dp) * dx
    end do
  end subroutine grid

  ! Subroutine to set initial conditions
  subroutine initial_cond()
    implicit none
    integer :: i
    call grid()

    allocate(u_global(0:nx), rho(0:nx), p(0:nx), c(0:nx))

    ! Left state (for x < 0.5)
    rho(0:nx/2) = 1.0d0
    p(0:nx/2) = 1.0d0
    u_global(0:nx/2) = 0.0d0

    ! Right state (for x >= 0.5)
    rho((nx/2 + 1):nx) = 0.1250d0
    p((nx/2 + 1):nx) = 0.1d0
    u_global((nx/2 + 1):nx) = 0.0d0

    ! Speed of sound calculation
    c = sqrt((gamma * p) / rho)

    ! Write initial conditions to file
    open(unit=10, file='shock_tube_output.dat', status='replace')
    write(10, '(A)') 'x, rho, p, u, c'
    do i = 0, nx
        write(10, '(F10.5, 3F10.5, F10.5)') x_center(i+1), rho(i), p(i), u_global(i), c(i)
    end do
    close(10)
  end subroutine initial_cond

  ! Subroutine to compute the time step
  subroutine time_step(c, dx, cfl, dt)
    real(dp), intent(in) :: c(:), dx, cfl
    real(dp), intent(out) :: dt
    dt = (cfl * dx) / maxval(c)
    print *, "Max wave speed:", maxval(c)
    print *, "Time step:", dt
  end subroutine time_step
 
  ! Subroutine to compute residuals
  subroutine residuals(F_left, F_right, dx, r)
    integer :: i
    real(dp), intent(in) :: F_left(3, nx), F_right(3, nx), dx
    real(dp), intent(out) :: r(3, nx)

    do i = 1, nx+1
      r(:,i) = (F_left(:,i) - F_right(:,i)) / dx
    end do
  end subroutine residuals

  ! Subroutine for time integration
  subroutine time_integration(U, r, dt)
    real(dp), intent(in out) :: U(3, nx+1), r(3, nx+1), dt
    integer :: i

    do i = 1, nx+1
      U(:,i) = U(:,i) - dt * r(:,i)
    end do
  end subroutine time_integration

  ! Subroutine to output the time step results
  subroutine output_time_step(U, x_center, time_step)
    real(dp), intent(in) :: U(3, nx+1), x_center(:)
    integer, intent(in) :: time_step
    integer :: i
    character(len=30) :: filename
    write(filename, '(A,I0,A)') 'time_step_', time_step, '.dat'
    open(unit=60, file=filename, status='replace')
    write(60, '(A)') 'x, rho, rho*u, E'
    do i = 0, nx
      write(60, '(F10.5, 3F10.5)') x_center(i+1), U(1,i+1), U(2,i+1), U(3,i+1)
    end do
    close(60)
  end subroutine output_time_step

end module calc_sub