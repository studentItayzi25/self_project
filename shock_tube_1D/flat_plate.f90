
! Start
!  |
!  V
! 1. Initialize Parameters and Domain
!    - Define grid points (x)
!    - Set initial conditions
!  |
!  V
! 2. Convert to Conserved Variables
!    - Calculate U_i from ρ, u, p
!  |
!  V
! 3. Compute Fluxes
!    - Calculate F(U)
!  |
!  V
! 4. Compute Numerical Flux
!    - Use Roe or HLL at interfaces
!  |
!  V
! 5. Update Conserved Variables
!    - Apply finite-volume scheme
!  |
!  V
! 6. Convert Back to Primitive Variables
!    - Calculate ρ, u, p from U
!  |
!  V
! 7. Check Termination
!    - t < T_final? Yes -> Repeat
!    - No -> Stop
!  |
!  V
! 8. Post-Processing
!    - Plot and compare results
!    - Compute error norms

module variables
implicit none

integer, parameter :: sp = kind(1.0)
integer, parameter :: dp = kind(1.0d0)
! grid requrment
integer, parameter :: nx = 1000



end module variables


module calc_sub
use variables

implicit none


contains

subroutine grid(x)
implicit none
real(dp) :: grid_points(nx + 1)
real(dp), parameter :: x_lbound = 0.0d0
real(dp), parameter :: x_rbound = 1.0d0
real(dp) :: dx
integer :: i 
character(len=20) :: filename 
integer :: x

dx = (x_rbound - x_lbound) / nx

write(filename, '(A,I0,A)') 'grid', nx, '.dat'
open(unit=10, file=filename, status='replace')
do i = 0, nx
grid_points(i+1) = x_lbound + float(i) * dx

write(10, *) grid_points(i+1)
end do

close(10)



 end subroutine grid

 end module calc_sub






program grid_test
  use calc_sub
  implicit none

  ! Declare variable x, which could be passed to grid
  integer :: x

  ! Call the grid subroutine to generate grid points
  call grid(x)

end program grid_test







! program flat_plate








! end program flat_plate