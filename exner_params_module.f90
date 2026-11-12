! =====================================================================
!  exner_params_module -- physical/numerical constants and run-time flags
!
!  Grid size and dx are not here: GeoClaw reads them at run time
!  (grid_module), so routines that need dx take it as an argument.
!  Run-time variables are overwritten from setprob.data by setprob.f90.
! =====================================================================
module exner_params_module

    implicit none
    save

    real(kind=8), parameter :: g   = 9.81d0
    real(kind=8), parameter :: pi  = 3.141592653589793d0
    real(kind=8), parameter :: eps = 1.0d-30          ! guards divisions

    ! ---- Channel geometry ----
    real(kind=8), parameter :: L_channel = 6.0d0     ! length [m]
    real(kind=8), parameter :: B_channel = 0.04d0    ! width  [m]

    ! ---- Sediment ----
    real(kind=8), parameter :: d50      = 0.0027d0   ! median grain diameter [m]
    real(kind=8), parameter :: rho_s    = 2550.0d0   ! sediment density [kg/m^3]
    real(kind=8), parameter :: rho_w    = 1000.0d0   ! water density    [kg/m^3]
    real(kind=8), parameter :: porosity = 0.4d0      ! bed porosity [-]
    real(kind=8), parameter :: s_rel    = (rho_s - rho_w) / rho_w

    real(kind=8), parameter :: mann_n = 0.025d0      ! Manning n [s/m^(1/3)]

    ! ---- Boundary and initial conditions ----
    real(kind=8), parameter :: Q_in       = 0.0002d0  ! inlet discharge [m^3/s]
    real(kind=8), parameter :: h_in       = 0.0116d0    ! inlet depth [m]
    real(kind=8), parameter :: qs_in      = 3.0d-5     ! inlet sediment flux per width [m^2/s]
    real(kind=8), parameter :: S0         = 0.044d0    ! rigid-floor slope [-]
    real(kind=8), parameter :: h_bed_init = 0.04d0     ! initial sediment thickness [m]
    real(kind=8), parameter :: h_bed_ds   = h_bed_init ! pinned downstream thickness [m]

    ! Initial sinusoidal bed perturbation
    logical,      parameter :: PERTURB_ON  = .true.
    real(kind=8), parameter :: A_perturb   = 4.0d-5    ! amplitude [m]
    real(kind=8), parameter :: lam_perturb = 0.5d0     ! wavelength [m]

    ! ---- Wet/dry thresholds and velocity cap ----
    real(kind=8), parameter :: h_wet     = 5.0d-4      ! cells with h <= h_wet are dry [m]
    real(kind=8), parameter :: c_dry     = 1.0d-3      ! wave speed at dry-dry interfaces [m/s]
    real(kind=8), parameter :: u_in_est  = Q_in / (B_channel * h_in)   ! inlet velocity [m/s]
    real(kind=8), parameter :: c_in_est  = sqrt(g * h_in)              ! inlet wave speed [m/s]
    real(kind=8), parameter :: u_max_cap = 3.0d0 * (u_in_est + 2.0d0 * c_in_est)   ! [m/s]

    ! ---- Run-time (setprob.data) ----
    real(kind=8) :: bed_alpha = 0.0d0     ! optional numerical bed diffusion coefficient [-]
    logical :: sediment_on = .true.       ! .false. freezes z and b (pure shallow water)
    logical :: friction_on = .true.       ! .false. disables Manning friction
    real(kind=8) :: mann_n_flow = mann_n  ! Manning n of the flow friction (sediment closure keeps mann_n)

contains

    ! Velocity from (h,hu): zero below h_wet, regularised near h_wet
    ! (avoids 0/0), and capped at u_max_cap.
    pure function compute_velocity(h_loc, hu_loc) result(u_loc)
        real(kind=8), intent(in) :: h_loc, hu_loc
        real(kind=8) :: u_loc
        if (h_loc <= h_wet) then
            u_loc = 0.0d0
        else
            u_loc = 2.0d0 * h_loc * hu_loc / (h_loc**2 + max(h_loc**2, h_wet**2))
            u_loc = sign(min(abs(u_loc), u_max_cap), u_loc)
        end if
    end function compute_velocity

end module exner_params_module
