module interpolations
    use variables
    use boundary_conditions
    implicit none

    contains

    subroutine compute_left_right_states
        implicit none
        integer :: i

        do i = 1, nx
            rho_left(i) = rho_c(i-1)
            u_left(i) = u_c(i-1)
            p_left(i) = p_c(i-1)
            c_left(i) = c_c(i-1)
            rho_right(i) = rho_c(i)
            u_right(i) = u_c(i)
            p_right(i) = p_c(i)
            c_right(i) = c_c(i)
        end do
    
    end subroutine compute_left_right_states

    subroutine U_left_right
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

    end subroutine U_left_right

end module interpolations

! module prim_con_in
!     use variables
!     use boundary_conditions
!     implicit none

!     contains

!     subroutine prim
!         implicit none
!         integer :: i
!         call boundary_conditions_init
!         do i = 1, nx
!             U(1, i) = rho_c(i)
!             U(2, i) = u_c(i)
!             U(3, i) = p_c(i)
!         end do
!     end subroutine prim

!     subroutine con
!         implicit none
!         integer :: i
!         call boundary_conditions_init
!         do i = 1, nx
!             F(1, i) = U(1, i) 
!             F(2, i) = U(2, i) * U(1, i)
!             ! E = p / (gamma - 1) + 0.5 * rho * u^2
!             F(3, i) = (U(3, i) / (gamma - 1.0d0)) + 0.5d0 * U(1, i) * U(2, i)**2 
!         end do
!     end subroutine con

!     subroutine prim_left_right
!         implicit none
!         integer :: i
!         call boundary_conditions_init
!         do i = 1, nx
!             U_l(1, i) = rho_left(i)
!             U_l(2, i) = u_left(i)
!             U_l(3, i) = p_left(i)
!             U_r(1, i) = rho_right(i)
!             U_r(2, i) = u_right(i)
!             U_r(3, i) = p_right(i)
!         end do
!     end subroutine prim_left_right

!     subroutine con_left_right
!         implicit none
!         integer :: i
!         call boundary_conditions_init
!         do i = 1, nx
!             F_l(1, i) = U_l(1, i) 
!             F_l(2, i) = U_l(2, i) * U_l(1, i)
!             F_l(3, i) = (U_l(3, i) / (gamma - 1.0d0)) + 0.5d0 * U_l(1, i) * U_l(2, i)**2 

!             F_r(1, i) = U_r(1, i) 
!             F_r(2, i) = U_r(2, i) * U_r(1, i)
!             F_r(3, i) = (U_r(3, i) / (gamma - 1.0d0)) + 0.5d0 * U_r(1, i) * U_r(2, i)**2 
!         end do
!     end subroutine con_left_right

! end module prim_con_in