module variables
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 500
  real(dp), parameter :: x_lbound = 0.0d0
  real(dp), parameter :: x_rbound = 1.0d0
  real(dp), parameter :: gamma = 1.4d0
  real(dp), parameter :: p_inf = 6.0d8
  real(dp), parameter :: cfl = 0.2d0
  real(dp), allocatable :: u_global(:), rho(:), p(:), c(:), x_center(:), x_edges(:)
  real(dp) :: dx, dt
end module variables

module calc_sub
  use variables
  implicit none
contains

  subroutine grid()
    implicit none
    integer :: i
    dx = (x_rbound - x_lbound) / nx
    ! Centers
    allocate(x_center(nx+1), x_edges(nx+1))
    do i = 0, nx-1
        x_center(i+1) = x_lbound + (real(i, dp) + 0.5d0) * dx
    end do
    ! Edges
    do i = 0, nx
        x_edges(i+1) = x_lbound + real(i, dp) * dx
    end do
  end subroutine grid
    allocate(u_global(nx+1), rho(nx+1), p(nx+1), c(nx+1))
    ! Initialize conditions based on scenario
    if (scenario == 1) then
        ! Example: Sod shock tube problem
        do i = 1, nx/2
            rho(i) = 1.0d0
            p(i) = 1.0d0
            u_global(i) = 0.0d0
        end do
        do i = nx/2+1, nx
            rho(i) = 0.125d0
            p(i) = 0.1d0
            u_global(i) = 0.0d0
        end do
    else if (scenario == 2) then
        ! Another scenario (e.g., different initial conditions)
        ! Initialize rho, p, u_global accordingly
    end if
    c = sqrt(gamma * (p + p_inf) / rho)
  end subroutine initial_cond

  subroutine compute_U_F(rho, vel, p, U, F)
    implicit none
    real(dp), intent(in) :: rho(:), vel(:), p(    do i = 1, nx
        E = p(i) / (gamma - 1.0) + 0.5 * rho(i) * vel(i)**2
        U(1, i) = rho(i)
        U(2, i) = rho(i) * vel(i)
        U(3, i) = E

        F(1, i) = rho(i) * vel(i)
        F(2, i) = rho(i) * vel(i)**2 + p(i)
        F(3, i) = vel(i) * (E + p(i))
    end do
  end subroutine compute_U_F

  subroutine compute_flux_hll(U, F_left, F_right, F_conv)
    use variables
    implicit none
    real(dp), intent(in) :: U(:,:), F_left(:,:), F_right(:,:)
    real(dp), intent(out) :: F_conv(:,:)
    integer :: i
    real(dp) :: rho_l, rho_r, u_l, u_r, p_l, p_r, c_l, c_r, E_l, E_r
    real(dp) :: S_L, S_R, S_star

    do i = 1, size(U, 2) - 1
        ! Extract left and right states
        rho_l = U(1, i)
        rho_r = U(1, i+1)
        u_l = U(2, i) / rho_l
        u_r = U(2, i+1) / rho_r
        E_l = U(3, i)
        E_r = U(3, i+1)
        p_l = (gamma - 1.0d0) * (E_l - 0.5d0 * rho_l * u_l**2)
        p_r = (gamma - 1.0d0) * (E_r - 0.5d0 * rho_r * u_r**2)
        c_l = sqrt(gamma * (p_l + p_inf) / rho_l)
        c_r = sqrt(gamma * (p_r + p_inf) / rho_r)

        ! Compute wave speeds
        S_L = min(u_l - c_l, u_r - c_r)
        S_R = max(u_l + c_l, u_r + c_r)

        if (S_L >= 0.0d0) then
            ! Entire flux from the left state
            F_conv(:, i) = F_left(:, i)
        else if (S_R <= 0.0d0) then
            ! Entire flux from the right state
            F_conv(:, i) = F_right(:, i)
        else
            ! HLL flux
            S_star = (p_r - p_l + rho_l * u_l * (S_L - u_l) - rho_r * u_r * (S_R - u_r)) / &
                     (rho_l * (S_L - u_l) - rho_r * (S_R - u_r))
            F_conv(:, i) = (S_R * F_left(:, i) - S_L * F_right(:, i) + S_L * S_R * &
                           (U(:, i+1) - U(:, i))) / (S_R - S_L)
        end if
    end do
  end subroutine compute_flux_hll

  subroutine residuals(F_left, F_right, dx, r)
    implicit none
    integer :: i
    real(dp), intent(in) :: F_left(3, nx+1), F_right(3, nx+1), dx
    real(dp), intent(out) :: r(3, nx+1)

    do i = 1, nx
        r(:, i) = (F_left(:, i) - F_right(:, i)) / dx
    end do
  end subroutine residuals

  subroutine time_integration(U, r, dt)
    implicit none
    real(dp), intent(in out) :: U(3, nx+1), r(3, nx+1), dt
    integer :: i

    do i = 1, nx
        U(:, i) = U(:, i) - dt * r(:, i)
    end do
  end subroutine time_integration

  subroutine output_time_step(U, x_center, time_step)
    implicit none
    real(dp), intent(in) :: U(3, nx+1), x_center(:)
    integer, intent(in) :: time_step
    integer :: i
    character(len=30) :: filename
    write(filename, '(A,I0,A)') 'time_step_', time_step, '.dat'
    open(unit=60, file=filename, status='replace')
    write(60, '(A)') 'x, rho, rho*u, E'
    do i = 1, nx
        write(60, '(F10.5, 3F10.5)') x_center(i), U(1,i), U(2,i), U(3,i)
    end do
    close(60)
  end subroutine output_time_step

end module calc_sub

program shock_tube
  use variables
  use calc_sub
  implicit none

  real(dp) :: U(3, nx+1), F_left(3, nx+1), F_right(3, nx+1), F_conv(3, nx+1), r(3, nx+1)
  integer :: scenario, time_step, max_steps
  real(dp) :: t, t_end, residual_normtime
  max_steps = 1000

  write(*,*) "Choose scenario: 1 (Eq. 5.2) or 2 (Eq. 5.3)"
  read(*,*) scenario

  call grid()
  call initial_cond(scenario)
  call compute_U_F(rho, u_global, p, U, F_left)

  time_step = 0
  do while (t < t_end .and. time_step < max_steps)
    ! Compute HLL fluxes
    call compute_flux_hll(U, F_left, F_right, F_conv)

    ! Compute residuals and update solution
    call residuals(F_conv, F_conv, dx, r)
    call time_integration(U, r, dt)

    ! Compute the magnitude of the residual
    residual_norm = sqrt(sum(reshape(r, (/size(r)/))**2))

    ! Print the time and the magnitude of the residual
    write(*, '(F10.5, 2X, E12.5)') t, residual_norm

    t = t + dt
    time_step = time_step + 1

    ! Output at specific intervals
    if (mod(time_step, 50) == 0) then
      call output_time_step(U, x_center, time_step)
    end if
  end do

  ! Final output
  call output_time_step(U, x_center, time_step)
end program shock_tube
