module grid
    use variables
    implicit none
    
    contains
    
    subroutine grid_init
        implicit none
        integer :: i

        ! Allocate arrays
        call allocate_variables
        dx = (x_r - x_l) / real(nx)
        
        do i = 0, nx
            x_face(i) = x_l + (i) * dx    
        end do

        do i = 1, nx 
            x_center(i) = 0.5d0 * (x_face(i) + x_face(i-1)) 
        end do

    end subroutine grid_init

end module grid