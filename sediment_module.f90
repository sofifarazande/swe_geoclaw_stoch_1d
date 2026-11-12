! =====================================================================
!  sediment_module -- stochastic bedload closure (Ancey): per-cell
!  equilibrium transport rate, particle velocity, entrainment rate and
!  activity diffusivity, plus the random numbers used by src1.
! =====================================================================
module sediment_module

    use exner_params_module
    implicit none
    save

    ! ---- Run-time (setprob.data) ----
    real(kind=8) :: cd_dep        = 0.06d0   ! sigma = cd_dep*ws/d50  (deposition rate)
    real(kind=8) :: mu_over_sigma = 0.85d0   ! mu = mu_over_sigma*sigma (entrainment rate)

    ! ---- Closure constants ----
    real(kind=8), parameter :: bagnold_a       = 4.2d0    ! particle-velocity coefficients
    real(kind=8), parameter :: bagnold_b_coeff = 0.11d0

    real(kind=8), parameter :: tau_star_0   = 0.035d0             ! critical Shields number, flat bed
    real(kind=8), parameter :: alpha_fric   = 47.0d0 * pi / 180.0d0   ! grain friction angle [rad]
    real(kind=8), parameter :: tau_star_max = 0.3d0               ! cap on the Shields number

    ! Guo/Recking equilibrium bedload law
    real(kind=8), parameter :: guo_c1 = 10.0d0
    real(kind=8), parameter :: guo_p1 = 16.0d0
    real(kind=8), parameter :: guo_c2 = 0.0768d0
    real(kind=8), parameter :: guo_p2 = 8.0d0
    real(kind=8), parameter :: guo_p3 = 13.0d0 / 8.0d0

    real(kind=8), parameter :: min_us    = 0.10d0       ! minimum particle velocity [m/s]
    real(kind=8), parameter :: h_sed_min = 5.0d0 * h_wet   ! no transport below this depth [m]
    real(kind=8), parameter :: sqrt_sgd  = sqrt(s_rel * g * d50)

    integer,      parameter :: n_noise_substeps = 25       ! Euler-Maruyama substeps per dt
    real(kind=8), parameter :: b_max_factor     = 20.0d0   ! b_max = b_max_factor*qs_eq_max/min_us
    real(kind=8), parameter :: noise_scale      = 1.0d0    ! noise amplitude factor (0 = no noise)

    ! Reaction block: fixed physical length (grain diameters), independent of dx
    real(kind=8), parameter :: n_diam_reaction = 15.0d0
    real(kind=8), parameter :: ell_reaction    = n_diam_reaction * d50
    real(kind=8), parameter :: v_p_grain       = (pi / 6.0d0) * d50**3   ! grain volume [m^3]

    ! ---- Set once by init_stochastic_transport ----
    real(kind=8) :: ws_g    = 0.0d0   ! settling velocity [m/s] (Dietrich 1982)
    real(kind=8) :: sigma_g = 0.0d0   ! deposition rate [1/s]
    real(kind=8) :: mu_g    = 0.0d0   ! collective entrainment rate [1/s]
    real(kind=8) :: kappa_g = 0.0d0   ! sigma - mu [1/s]
    real(kind=8) :: b_max_g = 0.0d0   ! ceiling on b [m]
    integer      :: block_ncells = 1  ! grid cells per reaction block

contains

    ! Must run after mu_over_sigma and cd_dep are read (see setprob.f90).
    subroutine init_stochastic_transport(dx_in)
        real(kind=8), intent(in) :: dx_in
        real(kind=8) :: Phi_max, qs_eq_max
        integer :: seed_size
        integer, allocatable :: seed_array(:)

        ! Fixed seed: reproducible runs, same noise across parameter sweeps.
        call random_seed(size=seed_size)
        allocate(seed_array(seed_size))
        seed_array = 20261007
        call random_seed(put=seed_array)

        ws_g    = ws_settle()
        sigma_g = cd_dep * ws_g / d50
        mu_g    = mu_over_sigma * sigma_g
        kappa_g = sigma_g - mu_g

        Phi_max   = (guo_c1 * tau_star_max)**guo_p1 * &
                    (1.0d0 + (tau_star_max / guo_c2)**guo_p2)**(-guo_p3)
        qs_eq_max = Phi_max * sqrt(s_rel * g * d50**3)
        b_max_g   = b_max_factor * qs_eq_max / min_us

        block_ncells = max(1, nint(ell_reaction / dx_in))
    end subroutine init_stochastic_transport

    ! Shields number from Manning friction, capped at tau_star_max.
    pure function shields(h_loc, u_loc) result(ts)
        real(kind=8), intent(in) :: h_loc, u_loc
        real(kind=8) :: ts, cf
        cf = g * mann_n**2 / (h_loc**(1.0d0/3.0d0) + eps)
        ts = cf * u_loc**2 / (s_rel * g * d50 + eps)
        ts = min(ts, tau_star_max)
    end function shields

    ! Equilibrium transport rate from the Shields number (Guo/Recking).
    pure function qs_guo(ts) result(qs_eq)
        real(kind=8), intent(in) :: ts
        real(kind=8) :: qs_eq, Phi_val
        Phi_val = (guo_c1 * ts)**guo_p1 * (1.0d0 + (ts / guo_c2)**guo_p2)**(-guo_p3)
        qs_eq   = Phi_val * sqrt(s_rel * g * d50**3)
    end function qs_guo

    ! Mean particle velocity (Bagnold-type), with a slope-dependent
    ! critical Shields number: slope_loc = -dz/dx > 0 (downhill) lowers
    ! tau_star_c (Fernandez Luque & van Beek 1976).
    pure function u_p_from_shields(ts, slope_loc) result(u_p)
        real(kind=8), intent(in) :: ts, slope_loc
        real(kind=8) :: u_p, tau_star_c
        tau_star_c = max(tau_star_0 * sin(alpha_fric - atan(slope_loc)) / sin(alpha_fric), 0.0d0)
        u_p = bagnold_b_coeff * ws_g + bagnold_a * sqrt_sgd * &
              (sqrt(max(ts, 0.0d0)) - sqrt(tau_star_c))
        u_p = max(u_p, min_us)
    end function u_p_from_shields

    ! Equilibrium transport rate qs_eq(h,u); zero in very shallow water.
    pure function qs_equilibrium(h_loc, u_loc) result(qs_eq)
        real(kind=8), intent(in) :: h_loc, u_loc
        real(kind=8) :: qs_eq
        if (h_loc <= h_sed_min) then
            qs_eq = 0.0d0
        else
            qs_eq = qs_guo(shields(h_loc, u_loc))
        end if
    end function qs_equilibrium

    ! Particle velocity u_p(h,u,slope); min_us in very shallow water.
    pure function particle_velocity(h_loc, u_loc, slope_loc) result(u_p)
        real(kind=8), intent(in) :: h_loc, u_loc, slope_loc
        real(kind=8) :: u_p
        if (h_loc <= h_sed_min) then
            u_p = min_us
        else
            u_p = u_p_from_shields(shields(h_loc, u_loc), slope_loc)
        end if
    end function particle_velocity

    ! Full per-cell closure: particle velocity u_p, entrainment rate
    ! lambda, and activity diffusivity Du.
    pure subroutine get_stochastic_parameters(h_loc, u_loc, slope_loc, &
                                              u_p, lambda_val, Du_val)
        real(kind=8), intent(in)  :: h_loc, u_loc, slope_loc
        real(kind=8), intent(out) :: u_p, lambda_val, Du_val
        real(kind=8) :: ts

        if (h_loc <= h_sed_min) then
            u_p        = min_us
            lambda_val = 0.0d0
            Du_val     = 0.0d0
            return
        end if

        ts         = shields(h_loc, u_loc)
        u_p        = u_p_from_shields(ts, slope_loc)
        lambda_val = max(qs_guo(ts) * kappa_g / u_p, 0.0d0)
        Du_val     = max(1.6d0 * abs(u_loc) / max(ws_g, eps) - 1.02d0, 0.05d0) * d50 * ws_g
    end subroutine get_stochastic_parameters

    ! Instantaneous transport rate qs = sign(u_p*b, u).
    pure function qs_stochastic(h_loc, u_loc, b_loc, u_p) result(qs_val)
        real(kind=8), intent(in) :: h_loc, u_loc, b_loc, u_p
        real(kind=8) :: qs_val
        if (h_loc <= h_wet .or. b_loc <= 0.0d0) then
            qs_val = 0.0d0
        else
            qs_val = sign(u_p * b_loc, u_loc)
        end if
    end function qs_stochastic

    ! Settling velocity, Dietrich (1982).
    pure function ws_settle() result(ws)
        real(kind=8) :: ws, nu_k, d_star, L, logF
        nu_k   = 1.0d-6
        d_star = s_rel * g * d50**3 / nu_k**2
        L      = log10(d_star)
        logF   = -3.76715d0 + 1.92944d0*L - 0.09815d0*L**2 - 0.00575d0*L**3 + 0.00056d0*L**4
        ws     = (s_rel * g * nu_k)**(1.0d0/3.0d0) * (10.0d0**logF)**(1.0d0/3.0d0)
    end function ws_settle

    ! Standard normal sample (Box-Muller).
    function normal_gaussian() result(gauss_val)
        real(kind=8) :: gauss_val, u1, u2
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, eps)
        gauss_val = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
    end function normal_gaussian

end module sediment_module
