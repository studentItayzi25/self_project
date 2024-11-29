module fluxes_methods
  use variables
  use calc_sub
  use ieee_arithmetic
  implicit none
contains

  ! Subroutine to compute left and right states
  subroutine compute_left_right_states(x_edges, U, F, F_left, F_right)
    real(dp), intent(in) :: x_edges(:), U(:,:), F(:,:)
    real(dp), intent(out) :: F_left(:,:), F_right(:,:)
    integer :: i
    real(dp) :: rho_l, rho_r, u_l, u_r, p_l, p_r, c_l, c_r, E_l, E_r

    do i = 1, size(x_edges) - 1
      ! Get left and right states from U
      rho_l = U(1, i)
      rho_r = U(1, i+1)
      u_l = U(2, i)
      u_r = U(2, i+1)
      E_l = U(3, i)
      E_r = U(3, i+1)

      ! Ensure that rho_l and rho_r are not zero to avoid division by zero
      ! if (rho_l <= 0.0d0 .or. rho_r <= 0.0d0) then
      !   print *, "Error: Non-positive density encountered at index ", i
      !   print *, "rho_l = ", rho_l, " rho_r = ", rho_r
      !   stop
      ! end if

      ! Calculate pressures
      p_l = (gamma - 1.0d0) * (E_l - 0.5d0 * rho_l * u_l**2)
      p_r = (gamma - 1.0d0) * (E_r - 0.5d0 * rho_r * u_r**2)

      ! ! Ensure that pressures are not negative
      ! if (p_l < 0.0d0 .or. p_r < 0.0d0) then
      !   print *, "Error: Negative pressure encountered at index ", i
      !   print *, "p_l = ", p_l, " p_r = ", p_r
      !   stop
      ! end if

      ! Calculate speed of sound
      c_l = sqrt(gamma * p_l / rho_l)
      c_r = sqrt(gamma * p_r / rho_r)

      ! Reflective boundary conditions for flux at the first and last index
      if (i == 1) then
        F_left(:, i) = F(:, i+1)   ! Reflective boundary (use right state for left flux)
        F_right(:, i) = F(:, i)    ! Use the left flux for the boundary right
      else if (i == size(x_edges) - 1) then
        F_left(:, i) = F(:, i)    ! Use the left flux for the boundary left
        F_right(:, i) = F(:, i-1) 
      else
        F_left(:, i) = 0.5d0 * (F(:, i) + F(:, i+1)) - 0.5d0 * (abs(c_l) + abs(u_r)) * (U(:, i+1) - U(:, i))
        F_right(:, i) = 0.5d0 * (F(:, i) + F(:, i+1)) + 0.5d0 * (abs(c_r) + abs(u_r)) * (U(:, i+1) - U(:, i))
      end if
    end do
  end subroutine compute_left_right_states



! Input: Left state (rho_L, u_L, p_L), Right state (rho_R, u_R, p_R)
! Output: Numerical flux F_GLF

! 1. Compute the left and right sound speeds:
!    c_L = sqrt(gamma * p_L / rho_L)
!    c_R = sqrt(gamma * p_R / rho_R)

! 2. Compute the maximum wave speed:
!    S = max(|u_L| + c_L, |u_R| + c_R)

! 3. Compute the flux vectors for left and right states:
!    F_L = [rho_L * u_L, rho_L * u_L^2 + p_L, u_L * (E_L + p_L)]
!    F_R = [rho_R * u_R, rho_R * u_R^2 + p_R, u_R * (E_R + p_R)]

! 4. Compute the GLF flux:
!    F_GLF = 0.5 * (F_L + F_R) - 0.5 * S * (U_R - U_L)

! Return: F_GLF

  ! Subroutine for GLF method

  subroutine GLF(x_edges, U, F, F_left, F_right, F_GLF)
    real(dp), intent(in) :: x_edges(:), U(:,:), F(:,:)
    real(dp), intent(out) :: F_left(:,:), F_right(:,:), F_GLF(:,:)
    integer :: i

    ! Compute left and right states
    call compute_left_right_states(x_edges, U, F, F_left, F_right)

    ! Compute the convective fluxes
    do i = 1, size(U, 2) - 1
        F_GLF(:, i) = 0.5d0 * (F_left(:, i) + F_right(:, i))
    end do
  end subroutine GLF



    ! Input: Left state (rho_L, u_L, p_L), Right state (rho_R, u_R, p_R)
    ! Output: Numerical flux F_HLL

    ! 1. Compute the left and right sound speeds:
    !    c_L = sqrt(gamma * p_L / rho_L)
    !    c_R = sqrt(gamma * p_R / rho_R)

    ! 2. Compute the wave speed bounds:
    !    S_L = min(u_L - c_L, u_R - c_R)
    !    S_R = max(u_L + c_L, u_R + c_R)

    ! 3. Compute the flux vectors for left and right states:
    !    F_L = [rho_L * u_L, rho_L * u_L^2 + p_L, u_L * (E_L + p_L)]
    !    F_R = [rho_R * u_R, rho_R * u_R^2 + p_R, u_R * (E_R + p_R)]

    ! 4. Compute the HLL flux:
    !    If S_L >= 0:
    !        F_HLL = F_L
    !    Else if S_R <= 0:
    !        F_HLL = F_R
    !    Else:
    !        F_HLL = (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L)

    ! Return: F_HLL

    ! Subroutine for HLL method
  subroutine HLL(x_edges, U, F, F_left, F_right, F_HLL, S_L, S_R)
    real(dp), intent(in) :: x_edges(:), U(:,:), F(:,:)
    real(dp), intent(out) :: F_left(:,:), F_right(:,:), F_HLL(:,:)
    real(dp), intent(out) :: S_L, S_R
    integer :: i
    real(dp) :: rho_l, rho_r, u_l, u_r, p_l, p_r, c_l, c_r, E_l, E_r

    ! Compute left and right states
    call compute_left_right_states(x_edges, U, F, F_left, F_right)

    ! Additional HLL flux calculation
    do i = 1, size(x_edges) - 1
        rho_l = U(1, i)
        rho_r = U(1, i+1)
        u_l = U(2, i)
        u_r = U(2, i+1)
        E_l = U(3, i)
        E_r = U(3, i+1)
        p_l = (gamma - 1.0d0) * (E_l - 0.5d0 * rho_l * u_l**2)
        p_r = (gamma - 1.0d0) * (E_r - 0.5d0 * rho_r * u_r**2)
        c_l = sqrt(gamma * p_l / rho_l)
        c_r = sqrt(gamma * p_r / rho_r)

        S_L = min(u_l - c_l, u_r - c_r)
        S_R = max(u_l + c_l, u_r + c_r)

        if (S_L >= 0.0) then
            F_HLL(:, i) = F_left(:, i)  ! Region 1: Entirely left state
        else if (S_R <= 0.0) then
            F_HLL(:, i) = F_right(:, i)  ! Region 3: Entirely right state
        else
            ! Region 2: Both waves contribute
            F_HLL(:, i) = (S_R * F_left(:, i) - S_L * F_right(:, i) + S_L * S_R * (U(:, i+1) - U(:, i))) / (S_R - S_L)
        end if
    end do
  end subroutine HLL

    ! Input: Left state (rho_L, u_L, p_L), Right state (rho_R, u_R, p_R)
    ! Output: Numerical flux F_HLLC

    ! 1. Compute the left and right sound speeds:
    ! c_L = sqrt(gamma * p_L / rho_L)
    ! c_R = sqrt(gamma * p_R / rho_R)

    ! 2. Compute the wave speed bounds:
    ! S_L = min(u_L - c_L, u_R - c_R)
    ! S_R = max(u_L + c_L, u_R + c_R)

    ! 3. Compute the contact wave speed:
    ! S_M = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
    !         (rho_L * (S_L - u_L) - rho_R * (S_R - u_R))

    ! 4. Compute the flux vectors for left and right states:
    ! F_L = [rho_L * u_L, rho_L * u_L^2 + p_L, u_L * (E_L + p_L)]
    ! F_R = [rho_R * u_R, rho_R * u_R^2 + p_R, u_R * (E_R + p_R)]

    ! 5. Compute the HLLC flux:
    ! If S_L >= 0:
    !     F_HLLC = F_L
    ! Else if S_R <= 0:
    !     F_HLLC = F_R
    ! Else if S_L < 0 < S_M:
    !     U_M = Compute intermediate state based on S_L, S_M, and U_L
    !     F_HLLC = F_L + S_L * (U_M - U_L)
    ! Else:
    !     U_M = Compute intermediate state based on S_R, S_M, and U_R
    !     F_HLLC = F_R + S_R * (U_M - U_R)

    ! Return: F_HLLC

  subroutine HLLC(x_edges, U, F, F_left, F_right, F_HLLC, S_L, S_R)
    real(dp), intent(in) :: x_edges(:), U(:,:), F(:,:)
    real(dp), intent(out) :: F_left(:,:), F_right(:,:), F_HLLC(:,:)
    real(dp), intent(out) :: S_L, S_R
    integer :: i
    real(dp) :: rho_l, rho_r, u_l, u_r, p_l, p_r, c_l, c_r, E_l, E_r
    real(dp) :: S_M, U_M(3)

    call compute_left_right_states(x_edges, U, F, F_left, F_right)

    do i = 1, size(x_edges) - 1
        rho_l = U(1, i)
        rho_r = U(1, i+1)
        u_l = U(2, i)
        u_r = U(2, i+1)
        E_l = U(3, i)
        E_r = U(3, i+1)
        p_l = (gamma - 1.0d0) * (E_l - 0.5d0 * rho_l * u_l**2)
        p_r = (gamma - 1.0d0) * (E_r - 0.5d0 * rho_r * u_r**2)
        c_l = sqrt(gamma * p_l / rho_l)
        c_r = sqrt(gamma * p_r / rho_r)

        S_L = min(u_l - c_l, u_r - c_r)
        S_R = max(u_l + c_l, u_r + c_r)

        S_M = (p_r - p_l + rho_l * u_l * (S_L - u_l) - rho_r * u_r * (S_R - u_r)) / &
            (rho_l * (S_L - u_l) - rho_r * (S_R - u_r))

        if (0.0d0 < S_M .and. S_M < S_L) then
            U_M(1) = rho_l * (S_L - u_l) / (S_L - S_M)
            U_M(3) = (S_M * (rho_l * E_l - rho_l * u_l * (S_L - u_l)) + p_l * (S_M - u_l)) / (S_L - S_M)
        else
            U_M(1) = rho_r * (S_R - u_r) / (S_R - S_M)
            U_M(3) = (S_M * (rho_r * E_r - rho_r * u_r * (S_R - u_r)) + p_r * (S_M - u_r)) / (S_R - S_M)
        end if
        U_M(2) = U_M(1) * S_M

        if (S_L >= 0.0) then
            F_HLLC(:, i) = F_left(:, i)
        else if (S_R <= 0.0) then
            F_HLLC(:, i) = F_right(:, i)
        else if (0.0d0 < S_L .and. S_M >= 0.0) then
            F_HLLC(:, i) = F_left(:, i)
        else
            F_HLLC(:, i) = F_HLLC(:, i) + S_L * (U_M - U(:, i))
        end if
    end do
  end subroutine HLLC

  



end module fluxes_methods
