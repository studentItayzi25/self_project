module boundary_conditions
    use variables
    implicit none

    contains

    subroutine boundary_conditions_init
        implicit none
  
        ! set boundary conditions for the left side
        rho_c(1:nx/2) = 1.0d0
        u_c(1:nx/2) = 0.0d0
        p_c(1:nx/2) = 1.0d0

        ! set boundary conditions for the right side
        rho_c((nx/2 + 1):nx) = 0.125d0
        u_c((nx/2 + 1):nx) = 0.0d0
        p_c((nx/2 + 1):nx) = 0.1d0

        c_c(1:nx) = sqrt(gamma * p_c(1:nx) / rho_c(1:nx))

    end subroutine boundary_conditions_init

end module boundary_conditions