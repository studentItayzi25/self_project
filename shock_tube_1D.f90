! structure for code

program shock_tube
  use variables
  use grid
  use boundary_conditions
  use prim_con_in
  use interpolations
  use fluxes
  use residuals
  use time
  
  implicit none

  integer :: i

  ! |------- pre-process----------|

  ! Initialize grid and allocate arrays
  call grid_init

  open(unit=20, file='grid.dat', status='replace')
  write(20, '(A)') 'x center, x face'
  do i = 0, nx
    write(20, '(F10.5, F10.5)') x_center(i), x_face(i)
  end do
  close(20)

  ! Initialize initial condition
  call boundary_conditions_init

  open(unit=30, file='initial_condition.dat', status='replace')
  write(30, '(A)') 'rho, u, p, c'
  do i = 0, nx
    write(30, '(F10.5, F10.5, F10.5, F10.5)') rho_c(i), u_c(i), p_c(i), c_c(i)
  end do
  close(30)


  ! |-------solver----------|

  ! Interpolation of left and right states
  call compute_left_right_states

  open(unit=40, file='left_right_states.dat', status='replace')
  write(40, '(A)') 'rho_left, u_left, p_left, c_left, rho_right, u_right, p_right, c_right'
  do i = 1, nx
    write(40, '(F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)') &
          rho_left(i), u_left(i), p_left(i), c_left(i), rho_right(i), u_right(i), p_right(i), c_right(i)
  end do
  close(40)

  call U_left_right

  open(unit=50, file='U.dat', status='replace')
  write(50, '(A)') 'rho, u, p'
  do i = 1, nx
    write(50, '(F10.5, F10.5, F10.5)') U(1, i), U(2, i), U(3, i)
  end do
  close(50)

  call con_left_right

  open(unit=60, file='F.dat', status='replace')
  write(60, '(A)') 'rho left, u left, E left, rho right, u right, E right'
  do i = 1, nx
    write(60, '(F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)') F_l(1, i), F_l(2, i), F_l(3, i), F_r(1, i), F_r(2, i), F_r(3, i)
  end do
  close(60)

  call simple_fluxes

  open(unit=70, file='F_flux.dat', status='replace')
  write(70, '(A)') 'rho flux, u flux, E flux'
  do i = 1, nx
    write(70, '(F10.5, F10.5, F10.5)') F_flux(1, i), F_flux(2, i), F_flux(3, i)
  end do
  close(70)
  
  call residuals_flux

open(unit=80, file='residuals.dat', status='replace')
write(80, '(A)') 'rho residual, u residual, E residual'
do i = 1, nx
    print *, 'Residuals:', r(1, i), r(2, i), r(3, i)  ! Debug print
    write(80, '(F10.5, F10.5, F10.5)') r(1, i), r(2, i), r(3, i)
end do
close(80)

  call time_advance

  open(unit=90, file='final_condition.dat', status='replace')
  write(90, '(A)') 'rho, u, p, c'
  do i = 1, nx
    write(90, '(F10.5, F10.5, F10.5, F10.5)') U_dt(1, i), U_dt(2, i), U_dt(3, i), c_c(i)
  end do
  close(90)


  

  ! Compute fluxes
  ! call GLF()
  ! call HLL()
  ! call HLLC()
  ! call residuals()
  ! call time_integration()

  ! |-------post-process----------|

end program shock_tube