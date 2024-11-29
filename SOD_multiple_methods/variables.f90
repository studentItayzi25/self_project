module variables
  implicit none
  ! Define precision and parameters
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 100
  real(dp), parameter :: x_lbound = 0.0d0
  real(dp), parameter :: x_rbound = 1.0d0
  real(dp), parameter :: gamma = 1.4d0
  real(dp), parameter :: cfl = 0.2d0
  real(dp), parameter :: Ma = 5.0d0
  integer :: method = 2  ! Define method as an integer variable
  ! Declare allocatable arrays
  real(dp), allocatable :: u_global(:), rho(:), p(:), c(:), x_center(:), x_edges(:)
  real(dp) :: dx, dt

contains

  ! Subroutine to print available methods and read user input
  subroutine select_method()
    print *, "Select a method:"
    print *, "1: GLF"
    print *, "2: HLL"
    print *, "3: HLLC"
    print *, "Enter the number corresponding to the desired method:"
    read(*, *) method
  end subroutine select_method

end module variables