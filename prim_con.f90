module prim_con_in
use variables
use boundary_conditions
implicit none

contains
subroutine prim
    implicit none
    integer :: i
    do i = 1, nx
        U(1, i) = rho_c(i)
        U(2, i) = rho_c(i) * u_c(i)
        U(3, i) = p_c(i) / (gamma - 1.0d0) + 0.5d0 * rho_c(i) * u_c(i)**2
    end do
end subroutine prim

subroutine con
    implicit none
    integer :: i
    do i = 1, nx
        F(1, i) = U(1, i)
        F(2, i) = U(2, i)
        F(3, i) = U(3, i)
    end do
end subroutine con

subroutine prim_left_right
    implicit none
    integer :: i
    do i = 1, nx
        U_l(1, i) = rho_left(i)
        U_l(2, i) = u_left(i)
        U_l(3, i) = p_left(i)
        U_r(1, i) = rho_right(i)
        U_r(2, i) = u_right(i)
        U_r(3, i) = p_right(i)
    end do
end subroutine prim_left_right

subroutine con_left_right
    implicit none
    integer :: i
    real(dp) :: E_left, E_right
    do i = 1, nx
        E_left = p_left(i) / (gamma - 1.0d0) + 0.5d0 * rho_left(i) * u_left(i)**2
        E_right = p_right(i) / (gamma - 1.0d0) + 0.5d0 * rho_right(i) * u_right(i)**2

        F_l(1, i) = rho_left(i)
        F_l(2, i) = rho_left(i) * u_left(i)
        F_l(3, i) = rho_left(i) * E_left

        F_r(1, i) = rho_right(i)
        F_r(2, i) = rho_right(i) * u_right(i)
        F_r(3, i) = rho_right(i) * E_right
    end do
end subroutine con_left_right

end module prim_con_in