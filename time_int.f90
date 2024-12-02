module time
use variables

contains

subroutine time_step
    implicit none
    real(dp) :: max_speed

    max_speed = maxval(abs(u_c(1:nx)) + c_c(1:nx))
    dt = cfl * dx / max_speed

    print *, 'Time step: ', dt
end subroutine time_step

subroutine time_advance
    implicit none
    integer :: i
    call time_step

    do i = 1, nx
        U_dt(1, i) = U(1, i) - dt * r(1, i)
        U_dt(2, i) = U(2, i) - dt * r(2, i)
        U_dt(3, i) = U(3, i) - dt * r(3, i)
    end do
end subroutine time_advance

end module time