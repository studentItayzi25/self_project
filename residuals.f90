module residuals
use variables
use fluxes
implicit none
contains

subroutine residuals_flux
    implicit none
    integer :: i
    call con_left_right

    do i = 1, nx
        r(1, i) =  (F_l(1, i) - F_r(1, i)) / dx
        r(2, i) =  (F_l(2, i) - F_r(2, i)) / dx
        r(3, i) =  (F_l(3, i) - F_r(3, i)) / dx
    end do

end subroutine residuals_flux

end module residuals