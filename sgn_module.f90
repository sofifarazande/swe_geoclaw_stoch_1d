! =====================================================================
!  sgn_module -- non-hydrostatic (Serre-Green-Naghdi) correction to hu,
!  adapted from GeoClaw's 1d_classic bouss_module (SGN branch, uniform
!  Cartesian grid, alpha = 1.153 improved dispersion).
!
!  Changes from GeoClaw: the bed is the evolving q(3) (not aux(1)); SGN
!  is active where the current depth exceeds sgn_min_depth on the whole
!  3-cell stencil; the correction is ramped smoothly to zero near the
!  outlet (bed pin); the bed in the ghost cells is extrapolated linearly.
!
!  Bed filter (sgn_bed_sigma > 0): the correction is computed for the
!  actual surface and discharge over the bed smoothed by a Gaussian of
!  standard deviation sgn_bed_sigma [m] (depth = eta - smoothed bed).
!  SGN is not valid for bed features shorter than the depth (measured
!  fixed beds with grain-scale roughness); the hydrostatic part (rp1,
!  src1) keeps the actual bed.
!
!  Each call solves (I + alpha T(h,B)) psi = rhs and applies
!     hu <- hu + dt * h * (g/alpha * eta_x - psi).
! =====================================================================
module sgn_module

    use exner_params_module, only: g
    implicit none

    ! ---- Run-time (setprob.data) ----
    logical      :: sgn_on        = .false.
    real(kind=8) :: sgn_bed_sigma = 0.0d0    ! bed filter std [m]; 0 = off

    real(kind=8), parameter :: alpha_sgn     = 1.153d0
    real(kind=8), parameter :: sgn_min_depth = 1.0d-3
    ! SGN weight 0 within sgn_off of the outlet (bed pin), rising (sin^2) to 1 at sgn_off + sgn_ramp
    real(kind=8), parameter :: sgn_off  = 0.05d0
    real(kind=8), parameter :: sgn_ramp = 0.10d0

contains

    subroutine sgn_correction(meqn, mbc, mx, dx, dt, q)
        integer, intent(in) :: meqn, mbc, mx
        real(kind=8), intent(in) :: dx, dt
        real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc)

        real(kind=8), dimension(1-mbc:mx+mbc) :: h, u, B, eta
        real(kind=8), dimension(0:mx+1) :: w, phi
        real(kind=8), dimension(mx) :: dl, dd, du, rhs
        real(kind=8) :: hr, Br, etar, Brr, wr, phir, ur, c2, c1, m, s
        real(kind=8) :: wgt(mx)
        logical :: active(mx)
        integer :: i

        do i = 1, mx
            s = (dble(mx) - dble(i) + 0.5d0) * dx
            wgt(i) = sin(0.5d0*acos(-1.0d0) * min(max((s - sgn_off)/sgn_ramp, 0.0d0), 1.0d0))**2
        end do
        c2 = 1.0d0 / dx**2
        c1 = 0.5d0 / dx

        B = q(3,:)
        call extrapolate_ghosts(B)
        h = q(1,:)
        eta = h + B
        if (sgn_bed_sigma > 0.0d0) then
            call smooth_bed(mbc, mx, dx, B)
            call extrapolate_ghosts(B)
            where (q(1,:) > sgn_min_depth) h = max(eta - B, sgn_min_depth)
        end if
        do i = 1-mbc, mx+mbc
            u(i) = 0.0d0
            if (h(i) > sgn_min_depth) u(i) = q(2,i) / h(i)
        end do

        do i = 0, mx+1
            Brr    = c2 * (B(i+1) - 2.0d0*B(i) + B(i-1))
            ur     = c1 * (u(i+1) - u(i-1))
            w(i)   = u(i)**2 * Brr
            phi(i) = ur**2
        end do

        do i = 1, mx
            active(i) = minval(h(i-1:i+1)) > sgn_min_depth .and. wgt(i) > 0.0d0
            if (.not. active(i)) then
                dl(i) = 0.0d0;  dd(i) = 1.0d0;  du(i) = 0.0d0;  rhs(i) = 0.0d0
                cycle
            end if
            hr   = c1 * (h(i+1) - h(i-1))
            Br   = c1 * (B(i+1) - B(i-1))
            etar = c1 * (eta(i+1) - eta(i-1))
            Brr  = c2 * (B(i+1) - 2.0d0*B(i) + B(i-1))
            wr   = c1 * (w(i+1) - w(i-1))
            phir = c1 * (phi(i+1) - phi(i-1))

            dd(i) = 1.0d0 + alpha_sgn * (2.0d0*h(i)**2/3.0d0 * c2 &
                                         + 0.5d0*h(i)*Brr + etar*Br)
            du(i) = alpha_sgn * (-h(i)**2/3.0d0 * c2 - h(i)*hr*c1)
            dl(i) = alpha_sgn * (-h(i)**2/3.0d0 * c2 + h(i)*hr*c1)
            rhs(i) = g/alpha_sgn * etar &
                     + 2.0d0*h(i) * (h(i)/3.0d0 * phir + phi(i)*(hr + 0.5d0*Br)) &
                     + 0.5d0*h(i)*wr + w(i)*etar
        end do

        ! Thomas algorithm; psi = 0 in the ghost cells
        do i = 2, mx
            m      = dl(i) / dd(i-1)
            dd(i)  = dd(i)  - m * du(i-1)
            rhs(i) = rhs(i) - m * rhs(i-1)
        end do
        rhs(mx) = rhs(mx) / dd(mx)
        do i = mx-1, 1, -1
            rhs(i) = (rhs(i) - du(i) * rhs(i+1)) / dd(i)
        end do

        do i = 1, mx
            if (active(i)) then
                etar = c1 * (eta(i+1) - eta(i-1))
                q(2,i) = q(2,i) + dt * wgt(i) * h(i) * (g/alpha_sgn * etar - rhs(i))
            end if
        end do
    contains

        subroutine extrapolate_ghosts(Bg)
            real(kind=8), intent(inout) :: Bg(1-mbc:mx+mbc)
            integer :: j
            do j = 0, 1-mbc, -1
                Bg(j) = 2.0d0*Bg(j+1) - Bg(j+2)
            end do
            do j = mx+1, mx+mbc
                Bg(j) = 2.0d0*Bg(j-1) - Bg(j-2)
            end do
        end subroutine extrapolate_ghosts

    end subroutine sgn_correction

    ! Gaussian smoothing of B(1:mx), kernel truncated at 3 sgn_bed_sigma.
    ! Beyond the ends B is continued linearly with the mean slope of the
    ! nearest kernel half-width, so the bed slope is not biased there.
    subroutine smooth_bed(mbc, mx, dx, B)
        integer, intent(in) :: mbc, mx
        real(kind=8), intent(in) :: dx
        real(kind=8), intent(inout) :: B(1-mbc:mx+mbc)

        real(kind=8), allocatable :: Be(:), wk(:)
        real(kind=8) :: sl, sr
        integer :: nk, i, n

        nk = max(1, ceiling(3.0d0 * sgn_bed_sigma / dx))
        n  = min(nk, mx - 1)
        allocate(Be(1-nk:mx+nk), wk(-nk:nk))
        Be(1:mx) = B(1:mx)
        sl = (B(1+n)  - B(1))      / dble(n)
        sr = (B(mx)   - B(mx-n))   / dble(n)
        do i = 1, nk
            Be(1-i)  = B(1)  - dble(i) * sl
            Be(mx+i) = B(mx) + dble(i) * sr
        end do
        do i = -nk, nk
            wk(i) = exp(-0.5d0 * (dble(i) * dx / sgn_bed_sigma)**2)
        end do
        wk = wk / sum(wk)
        do i = 1, mx
            B(i) = sum(wk * Be(i-nk:i+nk))
        end do
        deallocate(Be, wk)
    end subroutine smooth_bed

end module sgn_module
