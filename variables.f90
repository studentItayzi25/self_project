module variables
    implicit none
    
    ! Precision
    integer, parameter :: dp = kind(1.0d0)
    
    ! Integer parameters
    integer, parameter :: nx = 1000

    
    ! Real parameters
    real(dp), parameter :: x_l = 0.0d0, x_r = 1.0d0
    real(dp), parameter :: gamma = 1.4d0
    real(dp), parameter :: CFL = 0.2d0
    ! Real variables
    real(dp) :: dx
    real(dp) :: t, t_end
    real(dp) :: dt

    
    ! Allocatable arrays
    real(dp), allocatable :: x_center(:), x_face(:)
    real(dp), allocatable :: U(:,:), F(:,:), U_l(:,:), U_r(:,:), F_l(:,:), F_r(:,:), F_flux(:,:)
    real(dp), allocatable :: rho_c(:), u_c(:), p_c(:), c_c(:)
    real(dp), allocatable :: rho_left(:), u_left(:), p_left(:), c_left(:)
    real(dp), allocatable :: rho_right(:), u_right(:), p_right(:), c_right(:)
    real(dp), allocatable :: r(:,:)
    real(dp), allocatable :: U_dt(:,:)
contains
    subroutine allocate_variables
        implicit none
        allocate(x_center(1:nx), x_face(0:nx))
        allocate(U(1:3, 1:nx), F(1:3, 1:nx), U_l(1:3, 1:nx), U_r(1:3, 1:nx), F_l(1:3, 1:nx), F_r(1:3, 1:nx))
        allocate(rho_c(1:nx), u_c(1:nx), p_c(1:nx), c_c(1:nx))
        allocate(rho_left(1:nx), u_left(1:nx), p_left(1:nx), c_left(1:nx))
        allocate(rho_right(1:nx), u_right(1:nx), p_right(1:nx), c_right(1:nx))
        allocate(F_flux(1:3, 1:nx))
        allocate(r(1:3, 1:nx))
        allocate(U_dt(1:3, 1:nx))
    end subroutine allocate_variables
end module variables