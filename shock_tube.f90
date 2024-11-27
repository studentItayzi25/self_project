module variables
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 10
  real(dp), parameter :: x_lbound = 0.0d0
  real(dp), parameter :: x_rbound = 1.0d0
  real(dp), parameter :: gamma = 1.4d0
  real(dp), parameter :: cfl = 0.2d0
  real(dp), parameter :: Ma = 5.0d0
  real(dp), allocatable :: u_global(:), rho(:), p(:), c(:), x(:)
  real(dp) :: dx, dt
end module variables

module calc_sub
  use variables
  implicit none
contains

  subroutine grid()
    integer :: i
    dx = (x_rbound - x_lbound) / nx
    do i = 0, nx
      x(i+1) = x_lbound + real(i, dp) * dx
    end do
  end subroutine grid

  subroutine initial_cond()
    integer :: i
    call grid()

    allocate(u_global(0:nx), rho(0:nx), p(0:nx), c(0:nx))

    ! left state
    rho(0:nx/2) = 1.0d0
    p(0:nx/2) = 1.0d0
    u_global(0:nx/2) = 0.0d0

    ! right state
    rho((nx/2 + 1):nx) = 0.1250d0
    p((nx/2 + 1):nx) = 0.1d0
    u_global((nx/2 + 1):nx) = 0.0d0

    c = sqrt((gamma * p) / rho)

    open(unit=10, file='shock_tube_output.dat', status='replace')
    write(10, '(A)') 'x, rho, p, u, c'
    do i = 0, nx
      write(10, '(F10.5, 3F10.5, F10.5)') x(i+1), rho(i), p(i), u_global(i), c(i)
    end do
    close(10)
  end subroutine initial_cond

  subroutine time_step(c, dx, cfl, dt)
    real(dp), intent(in) :: c(:), dx, cfl
    real(dp), intent(out) :: dt
    dt = (cfl * dx) / maxval(c)
    print *, "Max wave speed:", maxval(c)
    print *, "Time step:", dt
  end subroutine time_step

  subroutine compute_U_F(rho, vel, p, U, F)
    ! Input variables
    real(dp), intent(in) :: rho(:), vel(:), p(:)
    ! Output variables
    real(dp), intent(out) :: U(:,:), F(:,:)
    ! Local variables
    real(dp) :: E
    integer :: i

    do i = 1, size(rho)
      ! Compute total energy per unit volume
      E = p(i) / (gamma - 1.0d0) + 0.5d0 * rho(i) * vel(i)**2

      ! Store conservative variables in U
      U(1, i) = rho(i)        ! rho
      U(2, i) = rho(i) * vel(i) ! rho u
      U(3, i) = E             ! Total energy

      ! Store fluxes in F
      F(1, i) = rho(i) * vel(i)               ! rho u
      F(2, i) = rho(i) * vel(i)**2 + p(i)     ! rho u^2 + p
      F(3, i) = vel(i) * (E + p(i))           ! u(E + p)
    end do
  end subroutine compute_U_F

end module calc_sub

program shock_tube
  use variables
  use calc_sub
  implicit none

  real(dp) :: U(3, nx+1), F(3, nx+1)
  integer :: i

  allocate(x(nx+1))

  ! Initialize the shock tube problem
  call initial_cond()

  ! Compute U and F
  call compute_U_F(rho, u_global, p, U, F)

  ! Write U and F to a file
  open(unit=20, file='conservative_and_flux.dat', status='replace')
  write(20, '(A)') 'x, rho, rho*u, E, rho*u, rho*u^2+p, u(E+p)'
  do i = 0, nx
    write(20, '(F10.5, 6F10.5)') x(i+1), U(1,i+1), U(2,i+1), U(3,i+1), F(1,i+1), F(2,i+1), F(3,i+1)
  end do
  close(20)

  ! Compute time step
  call time_step(c, dx, cfl, dt)

end program shock_tube
