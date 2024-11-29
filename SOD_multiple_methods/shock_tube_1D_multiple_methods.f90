program shock_tube
  use variables
  use fluxes_methods
  use ieee_arithmetic
  use calc_sub
  implicit none

  ! Declare variables
  real(dp) :: U(3, nx+1), F(3, nx+1), F_left(3, nx+1), F_right(3, nx+1), F_result(3, nx+1), r(3, nx+1)
  real(dp) :: S_L, S_R, max_wave_speed, local_dt
  integer :: i

  ! Allocate memory for arrays
  allocate(x_center(nx+1), x_edges(nx+1))

  ! Set initial conditions
  call initial_cond()

  ! Compute left and right states
  call compute_left_right_states(x_edges, U, F, F_left, F_right)

  ! Select the method based on user input
  call select_method()

  ! Call the appropriate flux computation method based on the method variable
  select case (method)
    case (1)
      call GLF(x_edges, U, F, F_left, F_right, F_result)
    case (2)
      call HLL(x_edges, U, F, F_left, F_right, F_result, S_L, S_R)
    case (3)
      call HLLC(x_edges, U, F, F_left, F_right, F_result, S_L, S_R)
    case default
      print *, 'Error: Unknown method ', method
      stop
  end select

  ! Compute residuals
  call residuals(F_left, F_right, dx, r)

  ! Perform time integration
  call time_integration(U, r, local_dt)

  ! Compute the next time step
  call time_step(c, dx, cfl, local_dt)

  ! Write conservative variables and fluxes to file
  open(unit=20, file='conservative_and_flux.dat', status='replace')
  write(20, '(A)') 'x, rho, rho*u, E, F1, F2, F3'
  do i = 0, nx
    write(20, '(F10.5, 6F10.5)') x_center(i+1), U(1,i+1), U(2,i+1), U(3,i+1), F(1,i+1), F(2,i+1), F(3,i+1)
  end do
  close(20)

  ! Write fluxes to file based on method
  select case (method)
    case (1)
      open(unit=30, file='fluxes_GLF.dat', status='replace')
    case (2)
      open(unit=30, file='fluxes_HLL.dat', status='replace')
    case (3)
      open(unit=30, file='fluxes_HLLC.dat', status='replace')
  end select
  write(30, '(A)') 'x, F_result, F_left, F_right'
  do i = 0, nx
    write(30, '(F10.5, 3F10.5)') x_edges(i+1), F_result(1,i+1), F_left(1,i+1), F_right(1,i+1)
  end do
  close(30)

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

  ! Calculate the maximum wave speed
  max_wave_speed = maxval(abs(u_global) + c)

  ! Print the maximum wave speed and time step only once
  print *, "Max wave speed: ", max_wave_speed
  print *, "Time step: ", local_dt

end program shock_tube