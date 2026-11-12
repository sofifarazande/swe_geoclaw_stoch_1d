! =====================================================================
!  src1 -- local source terms, Godunov-split after the flux step (rp1)
!
!    hu : -g h (dz/dx + Sf), bed slope + Manning friction (mann_n_flow),
!         then the SGN non-hydrostatic correction if sgn_on
!    z  : numerical bed diffusion (bed_alpha)
!    b  : diffusion Du b_xx, then entrainment/deposition drift +
!         birth-death noise, block-averaged over ell_reaction
!    z  : floor at aux(1) = z_base; last cell pinned to z_base + h_bed_ds
!
!  All updates read a snapshot of q taken on entry, so the explicit
!  update does not depend on loop order.
! =====================================================================
subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    use exner_params_module, only: g, h_wet, mann_n_flow, bed_alpha, h_bed_ds, &
                                    B_channel, u_in_est, c_in_est, compute_velocity, &
                                    sediment_on, friction_on
    use sediment_module, only: get_stochastic_parameters, kappa_g, mu_g, b_max_g, &
                                noise_scale, normal_gaussian, n_noise_substeps, &
                                block_ncells, v_p_grain
    use sgn_module, only: sgn_on, sgn_correction

    implicit none
    integer, intent(in) :: meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx,t,dt
    double precision, intent(in) :: aux(maux,1-mbc:mx+mbc)
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)

    double precision :: h_old(1-mbc:mx+mbc), hu_old(1-mbc:mx+mbc)
    double precision :: z_old(1-mbc:mx+mbc), b_old(1-mbc:mx+mbc)
    double precision :: u_a(mx), slope_a(mx), Du_a(mx), lambda_a(mx)
    double precision :: Sf, nu_z, u_p_dum
    double precision :: lambda_bar, k_factor, dt_sub, bm, bp, amp
    integer :: i, k, istart, iend, ncell_blk

    h_old  = q(1,:)
    hu_old = q(2,:)
    z_old  = q(3,:)
    b_old  = q(4,:)

    do i = 1, mx
        u_a(i)     = compute_velocity(h_old(i), hu_old(i))
        slope_a(i) = -(z_old(i+1) - z_old(i-1)) / (2.0d0*dx)   ! = -dz/dx
    end do

    ! ---- momentum: bed slope + friction ----
    do i = 1, mx
        if (h_old(i) > h_wet) then
            Sf = 0.0d0
            if (friction_on) Sf = mann_n_flow**2 * u_a(i) * abs(u_a(i)) / h_old(i)**(4.0d0/3.0d0)
            q(2,i) = hu_old(i) - dt * g * h_old(i) * (-slope_a(i) + Sf)
        end if
    end do

    if (sgn_on) call sgn_correction(meqn, mbc, mx, dx, dt, q)

    if (.not. sediment_on) return

    ! ---- bed diffusion ----
    nu_z = bed_alpha * dx * (u_in_est + c_in_est) / 4.0d0
    do i = 1, mx
        q(3,i) = z_old(i) + dt * nu_z * (z_old(i+1) - 2.0d0*z_old(i) + z_old(i-1)) / dx**2
    end do

    ! ---- b: diffusion ----
    do i = 1, mx
        call get_stochastic_parameters(max(h_old(i),h_wet), u_a(i), slope_a(i), &
                                       u_p_dum, lambda_a(i), Du_a(i))
        q(4,i) = b_old(i) + dt * Du_a(i) * (b_old(i+1) - 2.0d0*b_old(i) + b_old(i-1)) / dx**2
    end do

    ! ---- b: reaction + noise per block, full-truncation Euler-Maruyama
    ! substeps (b is floored only inside the coefficients, not in the state)
    if (dt > 0.0d0) then
        dt_sub = dt / dble(n_noise_substeps)
        do istart = 1, mx, block_ncells
            iend      = min(istart+block_ncells-1, mx)
            ncell_blk = iend - istart + 1

            lambda_bar = sum(lambda_a(istart:iend)) / dble(ncell_blk)
            k_factor   = v_p_grain / (dble(ncell_blk) * dx * B_channel)

            bm = sum(q(4,istart:iend)) / dble(ncell_blk)
            do k = 1, n_noise_substeps
                bp  = min(max(bm,0.0d0), b_max_g)
                amp = noise_scale * sqrt(max(2.0d0*mu_g*bp*k_factor, 0.0d0))
                bm  = bm + dt_sub*(lambda_bar - kappa_g*bp) + sqrt(dt_sub)*amp*normal_gaussian()
            end do
            q(4,istart:iend) = min(max(bm,0.0d0), b_max_g)
        end do
    end if

    ! ---- bed floor and downstream pin ----
    do i = 1, mx
        q(3,i) = max(q(3,i), aux(1,i))
    end do
    q(3,mx) = aux(1,mx) + h_bed_ds

end subroutine src1
