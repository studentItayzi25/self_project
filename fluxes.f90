module fluxes
use variables
use prim_con_in
implicit none

contains
subroutine simple_fluxes
    implicit none
    integer :: i
    call con_left_right

    do i = 1, nx
        F_flux(1, i) =  (F_l(1, i) - F_r(1, i)) / 2.0d0
        F_flux(2, i) =  (F_l(2, i) - F_r(2, i)) / 2.0d0
        F_flux(3, i) =  (F_l(3, i) - F_r(3, i)) / 2.0d0
    end do

end subroutine simple_fluxes

end module fluxes
