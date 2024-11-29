module variables
  implicit none
  ! Define precision and parameters
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 300
  real(dp), parameter :: x_lbound = 0.0d0
  real(dp), parameter :: x_rbound = 1.0d0
  real(dp), parameter :: gamma = 1.4d0
  real(dp), parameter :: cfl = 0.2d0
  real(dp), parameter :: Ma = 5.0d0
  integer, parameter :: method = 2
  ! Declare allocatable arrays
  real(dp), allocatable :: u_global(:), rho(:), p(:), c(:), x_center(:), x_edges(:)
  real(dp) :: dx, dt
end module variables

module calc_sub
  use variables
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

  ! Subroutine to compute conservative variables and fluxes
  subroutine compute_U_F(rho, vel, p, U, F, F_left, F_right, F_conv)
    real(dp), intent(in) :: rho(:), vel(:), p(:)
    real(dp), intent(out) :: U(:,:), F(:,:), F_conv(:,:), F_left(:,:), F_right(:,:)
    real(dp) :: E
    integer :: i

    do i = 1, size(rho)
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

    ! Compute the convective fluxes
    do i = 1, size(U, 2) - 1
        F_conv(:, i) = 0.5d0 * (F_left(:, i) + F_right(:, i))
    end do
  end subroutine compute_U_F

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

module fluxes_methods
  use variables
  implicit none
contains

  ! Subroutine for flux vector splitting
  subroutine flux_vector_splitting(U, F, F_left, F_right, F_conv, F_HLL, S_L, S_R)
    implicit none
    real(dp), intent(in) :: U(:,:), F(:,:)
    real(dp), intent(out) :: F_left(:,:), F_right(:,:), F_conv(:,:), F_HLL(:,:)
    real(dp), intent(out) :: S_L, S_R
    integer :: i
    real(dp) :: rho_l, rho_r, u_l, u_r, p_l, p_r, c_l, c_r, E_l, E_r

    do i = 1, size(U, 2) - 1
        ! Get left and right states from U
        rho_l = U(1, i)
        rho_r = U(1, i+1)
        u_l = U(2, i)
        u_r = U(2, i+1)
        E_l = U(3, i)
        E_r = U(3, i+1)
        p_l = (gamma - 1.0d0) * (E_l - 0.5d0 * rho_l * u_l**2)
        p_r = (gamma - 1.0d0) * (E_r - 0.5d0 * rho_r * u_r**2)
        c_l = sqrt(gamma * p_l / rho_l)
        c_r = sqrt(gamma * p_r / rho_r)

        ! Reflective boundary conditions for flux at the first and last index
        if (i == 1) then
            F_left(:, i) = F(:, i+1)   ! Reflective boundary (use right state for left flux)
            F_right(:, i) = F(:, i)    ! Use the left flux for the boundary right
        else if (i == size(U, 2) - 1) then
            F_left(:, i) = F(:, i)    ! Use the left flux for the boundary left
            F_right(:, i) = F(:, i-1) 
        else
            F_left(:, i) = 0.5d0 * (F(:, i) + F(:, i+1)) - 0.5d0 * (abs(c_l) + abs(u_r)) * (U(:, i+1) - U(:, i))
            F_right(:, i) = 0.5d0 * (F(:, i) + F(:, i+1)) + 0.5d0 * (abs(c_r) + abs(u_r)) * (U(:, i+1) - U(:, i))
        end if
    end do

    ! Additional HLL flux calculation
    do i = 1, size(U, 2) - 1
        S_L = min(u_l - c_l, u_r - c_r)
        S_R = max(u_l + c_l, u_r + c_r)

        if (S_L >= 0.0) then
            F_HLL(:, i) = F_left(:, i)  ! Region 1: Entirely left state
        else if (S_R <= 0.0) then
            F_HLL(:, i) = F_right(:, i)  ! Region 3: Entirely right state
        else
            ! Region 2: Both waves contribute
            F_HLL(:, i) = (S_R * F_left(:, i) - S_L * F_right(:, i) + S_L * S_R * (U(:, i+1) - U(:, i))) / (S_R - S_L)
        end if
    end do
  end subroutine flux_vector_splitting

end module fluxes_methods

program shock_tube
  use variables
  use calc_sub
  use fluxes_methods
  implicit none

  ! Declare variables
  real(dp) :: U(3, nx+1), F(3, nx+1), F_left(3, nx+1), F_right(3, nx+1), F_conv(3, nx+1), r(3, nx+1), F_HLL(3, nx+1)
  real(dp) :: S_L, S_R
  integer :: i

  ! Allocate memory for arrays
  allocate(x_center(nx+1), x_edges(nx+1))

  ! Set initial conditions
  call initial_cond()

  ! Compute conservative variables and fluxes
  call compute_U_F(rho, u_global, p, U, F, F_left, F_right, F_conv)

  ! Perform flux vector splitting
  call flux_vector_splitting(U, F, F_left, F_right, F_conv, F_HLL, S_L, S_R)

  ! Compute residuals
  call residuals(F_left, F_right, dx, r)

  ! Perform time integration
  call time_integration(U, r, dt)

  ! Output the time step results
  call output_time_step(U, x_center, 0)

  ! Compute the next time step
  call time_step(c, dx, cfl, dt)

  ! Write conservative variables and fluxes to file
  open(unit=20, file='conservative_and_flux.dat', status='replace')
  write(20, '(A)') 'x, rho, rho*u, E, rho*u, rho*u^2+p, u(E+p)'
  do i = 0, nx
    write(20, '(F10.5, 6F10.5)') x_center(i+1), U(1,i+1), U(2,i+1), U(3,i+1), F(1,i+1), F(2,i+1), F(3,i+1)
  end do
  close(20)

  ! Write fluxes to file based on method
  if (method == 1) then
    open(unit=30, file='fluxes.dat', status='replace')
    write(30, '(A)') 'x, F_conv, F_left, F_right'
    do i = 0, nx
      write(30, '(F10.5, 3F10.5)') x_edges(i+1), F_conv(1,i+1), F_left(1,i+1), F_right(1,i+1)
    end do
    close(30)
  else
    open(unit=30, file='fluxes_HLL.dat', status='replace')
    write(30, '(A)') 'x, F_hll, F_left, F_right'
    do i = 0, nx
      write(30, '(F10.5, 3F10.5)') x_edges(i+1), F_HLL(1,i+1), F_left(1,i+1), F_right(1,i+1)
    end do
    close(30)
  end if

  ! Write residuals to file
  open(unit=40, file='residuals.dat', status='replace')
  write(40, '(A)') 'x, residual_rho, residual_rho*u, residual_E'
  do i = 0, nx
    write(40, '(F10.5, 3F10.5)') x_center(i+1), r(1,i+1), r(2,i+1), r(3,i+1)
  end do
  close(40)

  ! Write final results to file
  open(unit=50, file='final.dat', status='replace')
  write(50, '(A)') 'x, rho, rho*u, E'
  do i = 0, nx
    write(50, '(F10.5, 3F10.5)') x_center(i+1), U(1,i+1), U(2,i+1), U(3,i+1)
  end do
  close(50)

  ! Compute the next time step
  call time_step(c, dx, cfl, dt)
end program shock_tube