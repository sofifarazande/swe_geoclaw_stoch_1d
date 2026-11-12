#!/usr/bin/env python3
"""
linear_stability_4eq.py - linear growth rate and migration speed versus
wavelength for the 4-equation (h,u,z,b) mean-field system, in which b
relaxes towards equilibrium at the finite rate kappa (linear_stability.py
instead uses qs = qs_eq(h,u) instantly).

Model (as solved by rp1/src1):
  b_t + (u_p*b)_x = lambda(h,u,S) - kappa*b + Du*b_xx    (mean field, no noise)
  (1-poro) z_t + (u_p*b)_x = 0
  u_p(h,u,S), lambda(h,u,S) respond instantly; S = -z_x is the local slope
  (base slope S0); b0 = lambda0/kappa, so that u_p0*b0 = qs_eq(h0,u0)

Linear system for (h',u',z',b') ~ exp(i k x + sigma t), with
qs' = u_p0 b' + b0 (Up_h h' + Up_u u' - ik Up_S z'):

  row h': [-ik u0,                -ik h0,               0,                               0                       ]
  row u': [-ik g + g*4 S0/(3 h0), -ik u0 - g*2 S0/u0,   -ik g,                           0                       ]
  row z': [-ik b0 Up_h/(1-poro),  -ik b0 Up_u/(1-poro), -(D_bed + b0 Up_S/(1-poro)) k^2, -ik up0/(1-poro)        ]
  row b': [Lam_h - ik b0 Up_h,    Lam_u - ik b0 Up_u,   -ik Lam_S - b0 Up_S k^2,         -ik up0 - kappa - Du k^2]

--sgn replaces the u' row by the linearised Serre (SGN) momentum equation
  D u + g eta_x + friction = (h0^2/3) D u_xx - (h0/2) D^2 z_x,   D = d/dt + u0 d/dx,
whose sigma^2 z term is eliminated with the z' row.

Branches: two gravity waves, one bed/Exner-like and one activity-like
(decaying at about kappa), told apart by their eigenvectors.

All constants are read from the solver (model_params.py); mu_over_sigma,
cd_dep, bed_alpha and dx default to the values in setprob.data/claw.data.
The base flow is the normal flow for Q_in unless --from-output (mean
late-time flow of a run) or --h0/--u0 are given.

Usage:
    python scripts/linear_stability_4eq.py [--h0 H --u0 U | --from-output DIR]
        [--mu-over-sigma M] [--cd-dep C] [--bed-alpha A] [--sgn] [--save FILE.png]
    python scripts/linear_stability_4eq.py --sweep-kappa [--save FILE.png]
    python scripts/linear_stability_4eq.py --sweep-bed-alpha [--save FILE.png]
"""
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from model_params import load, base_flow, settling_velocity

P = load()
g, mann_n, d50, s_rel = P['g'], P['mann_n'], P['d50'], P['s_rel']
porosity, S0, h_sed_min, eps = P['porosity'], P['S0'], P['h_sed_min'], P['eps']
tau_star_max, tau_star_0, alpha_fric = P['tau_star_max'], P['tau_star_0'], P['alpha_fric']
bagnold_a, bagnold_b_coeff, min_us = P['bagnold_a'], P['bagnold_b_coeff'], P['min_us']
guo_c1, guo_p1, guo_c2, guo_p2, guo_p3 = (P[k] for k in ('guo_c1', 'guo_p1', 'guo_c2',
                                                          'guo_p2', 'guo_p3'))
WS_G = settling_velocity(P)


def kappa_of(mu_over_sigma, cd_dep):
    """(kappa, sigma, mu) as set by init_stochastic_transport."""
    sigma_g = cd_dep * WS_G / d50
    mu_g = mu_over_sigma * sigma_g
    return sigma_g - mu_g, sigma_g, mu_g


def qs_equilibrium(h, u):
    """As sediment_module.f90's qs_equilibrium."""
    if h <= h_sed_min:
        return 0.0
    cf = g * mann_n**2 / (abs(h)**(1.0/3.0) + eps)
    tau_star = min(cf * u**2 / (s_rel*g*d50 + eps), tau_star_max)
    Phi = (guo_c1*tau_star)**guo_p1 * (1.0 + (tau_star/guo_c2)**guo_p2)**(-guo_p3)
    return Phi * np.sqrt(s_rel*g*d50**3)


def particle_velocity(h, u, slope=S0):
    """Particle velocity u_p, as sediment_module.f90's particle_velocity."""
    if h <= h_sed_min:
        return min_us
    sqrt_sgd = np.sqrt(s_rel*g*d50)
    cf = g * mann_n**2 / (h**(1.0/3.0) + eps)
    tau_star = min(cf * u**2 / (s_rel*g*d50 + eps), tau_star_max)
    theta = np.arctan(slope)
    tau_star_c = max(tau_star_0 * np.sin(alpha_fric-theta) / np.sin(alpha_fric), 0.0)
    u_p = bagnold_b_coeff*WS_G + bagnold_a*sqrt_sgd*(np.sqrt(max(tau_star, 0.0)) - np.sqrt(tau_star_c))
    return max(u_p, min_us)


def lambda_entrainment(h, u, kappa, slope=S0):
    """Entrainment rate lambda, as get_stochastic_parameters."""
    if h <= h_sed_min:
        return 0.0
    u_p = particle_velocity(h, u, slope)
    qs_eq = qs_equilibrium(h, u)
    return max(qs_eq * kappa / u_p, 0.0)


def activity_diffusivity(h, u):
    """Activity diffusivity Du, as get_stochastic_parameters."""
    if h <= h_sed_min:
        return 0.0
    return max(1.6*abs(u)/max(WS_G, eps) - 1.02, 0.05) * d50 * WS_G


def base_state(h0, u0, kappa):
    """u_p0, lambda0, b0 = lambda0/kappa and Du0 at (h0,u0)."""
    u_p0 = particle_velocity(h0, u0)
    lambda0 = lambda_entrainment(h0, u0, kappa)
    b0 = lambda0 / kappa
    Du0 = activity_diffusivity(h0, u0)
    return u_p0, lambda0, b0, Du0


def derivatives(h0, u0, kappa):
    """Central differences of u_p and lambda with respect to h, u and the
    slope S (at S0)."""
    eps_u = max(1.0e-4*max(abs(u0), 0.01), 1.0e-6)
    eps_h = max(1.0e-4*h0, 1.0e-8)
    eps_s = 1.0e-6

    Up_u = (particle_velocity(h0, u0+eps_u) - particle_velocity(h0, u0-eps_u)) / (2*eps_u)
    Up_h = (particle_velocity(h0+eps_h, u0) - particle_velocity(h0-eps_h, u0)) / (2*eps_h)
    Up_S = (particle_velocity(h0, u0, S0+eps_s) - particle_velocity(h0, u0, S0-eps_s)) / (2*eps_s)

    Lam_u = (lambda_entrainment(h0, u0+eps_u, kappa) - lambda_entrainment(h0, u0-eps_u, kappa)) / (2*eps_u)
    Lam_h = (lambda_entrainment(h0+eps_h, u0, kappa) - lambda_entrainment(h0-eps_h, u0, kappa)) / (2*eps_h)
    Lam_S = (lambda_entrainment(h0, u0, kappa, S0+eps_s)
             - lambda_entrainment(h0, u0, kappa, S0-eps_s)) / (2*eps_s)

    return Up_h, Up_u, Up_S, Lam_h, Lam_u, Lam_S


def L_matrix(k, h0, u0, kappa, u_p0, b0, Up_h, Up_u, Up_S, Lam_h, Lam_u, Lam_S,
             D_bed=0.0, Du0=0.0, sgn=False):
    """4x4 linear operator. D_bed: bed diffusivity [m^2/s] (src1.f90's
    nu_z); Du0: activity diffusivity [m^2/s]; sgn: Serre momentum row."""
    ik = 1j*k
    M = np.array([
        [-ik*u0,                    -ik*h0,                     0.0,                                          0.0],
        [-ik*g + g*4*S0/(3*h0),     -ik*u0 - g*2*S0/u0,         -ik*g,                                        0.0],
        [-ik*b0*Up_h/(1-porosity),  -ik*b0*Up_u/(1-porosity),   -(D_bed + b0*Up_S/(1-porosity))*k**2,         -ik*u_p0/(1-porosity)],
        [Lam_h - ik*b0*Up_h,        Lam_u - ik*b0*Up_u,         -ik*Lam_S - b0*Up_S*k**2,                     -ik*u_p0 - kappa - Du0*k**2],
    ], dtype=complex)
    if sgn:
        # (1+K)(sigma + ik u0) u' + (h0/2) ik (sigma + ik u0)^2 z' = hydrostatic terms,
        # with sigma z' = M[2].X and sigma^2 z' = sum_j M[2,j] sigma X_j
        K = (h0*k)**2 / 3.0
        e_u, e_z = np.eye(4)[1], np.eye(4)[2]
        rhs = M[1] + ik*u0*e_u - (1+K)*ik*u0*e_u
        s2 = sum(M[2, j]*M[j] for j in (0, 2, 3))
        rhs = rhs - 0.5*h0*ik*(s2 + 2*ik*u0*M[2] - (k*u0)**2*e_z)
        M[1] = rhs / ((1+K) + 0.5*h0*ik*M[2, 1])
    return M


def bed_diffusivity_from_bed_alpha(bed_alpha, dx):
    """nu_z = bed_alpha*dx*(u_in_est+c_in_est)/4, as in src1.f90."""
    return bed_alpha * dx * (P['u_in_est'] + P['c_in_est']) / 4.0


def bed_alpha_from_bed_diffusivity(D_bed, dx):
    return 4.0 * D_bed / (dx * (P['u_in_est'] + P['c_in_est']))


def dispersion_curve(h0, u0, kappa, wavelengths, D_bed=0.0, sgn=False):
    """Eigenvalues sorted by migration speed at each k. With the
    eigenvectors scaled by (h0, u0, h0, b0), the bed branch is the mode with
    the largest z' share and the activity branch, of the others, the one
    with the largest b' share."""
    u_p0, lambda0, b0, Du0 = base_state(h0, u0, kappa)
    Up_h, Up_u, Up_S, Lam_h, Lam_u, Lam_S = derivatives(h0, u0, kappa)
    scale = np.array([h0, abs(u0), h0, b0])[:, None]

    n = len(wavelengths)
    growth = np.zeros((n, 4))
    speed = np.zeros((n, 4))
    bed_branch = np.zeros(n, dtype=int)
    act_branch = np.zeros(n, dtype=int)
    for i, wl in enumerate(wavelengths):
        k = 2*np.pi / wl
        M = L_matrix(k, h0, u0, kappa, u_p0, b0, Up_h, Up_u, Up_S, Lam_h, Lam_u, Lam_S,
                     D_bed, Du0, sgn)
        eig4, V = np.linalg.eig(M)
        order = np.argsort(-eig4.imag / k)
        growth[i] = eig4.real[order]
        speed[i] = -eig4.imag[order] / k
        W = np.abs(V[:, order] / scale)
        W /= W.sum(axis=0)
        bed_branch[i] = np.argmax(W[2])
        others = [j for j in range(4) if j != bed_branch[i]]
        act_branch[i] = others[np.argmax(W[3, others])]

    return growth, speed, bed_branch, act_branch, u_p0, lambda0, b0, Du0, Up_h, Up_u, Lam_h, Lam_u


def run_one(h0, u0, mu_over_sigma, cd_dep, wavelengths, verbose=True, D_bed=0.0, sgn=False):
    kappa, sigma_g, mu_g = kappa_of(mu_over_sigma, cd_dep)
    growth, speed, bed_branch, act_branch, u_p0, lambda0, b0, Du0, Up_h, Up_u, Lam_h, Lam_u = \
        dispersion_curve(h0, u0, kappa, wavelengths, D_bed, sgn)
    idx = np.arange(len(wavelengths))
    bed_growth = growth[idx, bed_branch]
    bed_speed = speed[idx, bed_branch]
    act_growth = growth[idx, act_branch]

    if verbose:
        print(f"\nmu_over_sigma={mu_over_sigma:.4f}  cd_dep={cd_dep:.4f}  D_bed={D_bed:.4e} m^2/s  "
              f"Du0={Du0:.4e} m^2/s  "
              f"-> sigma_g={sigma_g:.4f} 1/s, mu_g={mu_g:.4f} 1/s, kappa={kappa:.5f} 1/s "
              f"(relaxation time 1/kappa={1/kappa:.3f} s)")
        print(f"  Base state: u_p0={u_p0*1000:.2f} mm/s, lambda0={lambda0:.4e}, "
              f"b0={b0:.4e} m  (check: u_p0*b0={u_p0*b0:.4e} vs qs_eq(h0,u0)="
              f"{qs_equilibrium(h0,u0):.4e})")
        print(f"  Up_h={Up_h:.4e}, Up_u={Up_u:.4e}, Lam_h={Lam_h:.4e}, Lam_u={Lam_u:.4e}")
        i_peak = np.argmax(bed_growth)
        if np.max(bed_growth) > 0:
            print(f"  Bed/Exner branch: fastest-growing wavelength "
                  f"{wavelengths[i_peak]*100:.2f} cm (growth {bed_growth[i_peak]:.4f} 1/s, "
                  f"speed {bed_speed[i_peak]*1000:.3f} mm/s "
                  f"{'UPSTREAM' if bed_speed[i_peak]<0 else 'downstream'})")
        else:
            print(f"  Bed/Exner branch: NO growing wavelength - least-decaying "
                  f"{bed_growth[i_peak]:.4e} 1/s at {wavelengths[i_peak]*100:.2f} cm")
        print(f"  Activity branch growth range: [{act_growth.min():.3f}, {act_growth.max():.3f}] 1/s "
              f"(expect ~ -kappa = {-kappa:.3f} since it mostly just relaxes toward equilibrium)")

    return dict(kappa=kappa, growth=growth, speed=speed, bed_branch=bed_branch,
                act_branch=act_branch, bed_growth=bed_growth, bed_speed=bed_speed,
                act_growth=act_growth)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--h0', type=float, help='base flow depth (m)')
    parser.add_argument('--u0', type=float, help='base flow velocity (m/s)')
    parser.add_argument('--from-output', metavar='DIR',
                         help='use the mean flow of the 2nd half of this run as base flow')
    parser.add_argument('--mu-over-sigma', type=float, default=P['mu_over_sigma'])
    parser.add_argument('--cd-dep', type=float, default=P['cd_dep'])
    parser.add_argument('--sgn', action='store_true',
                         help='Serre (SGN) momentum instead of hydrostatic shallow water')
    parser.add_argument('--bed-alpha', type=float, default=P['bed_alpha'],
                         help='numerical bed diffusion, converted to D_bed with src1.f90 '
                              'nu_z formula at --dx (default: setprob.data)')
    parser.add_argument('--dx', type=float, default=P['dx'],
                         help='grid spacing, used for D_bed and Du0 (default: claw.data)')
    parser.add_argument('--wl-min', type=float, default=0.02)
    parser.add_argument('--wl-max', type=float, default=0.60)
    parser.add_argument('--n', type=int, default=400)
    parser.add_argument('--sweep-kappa', action='store_true',
                         help='overlay several mu_over_sigma values (0.01..0.995) at fixed cd_dep')
    parser.add_argument('--sweep-bed-alpha', action='store_true',
                         help='overlay several --bed-alpha values at fixed mu_over_sigma and cd_dep')
    parser.add_argument('--save', metavar='FILE.png')
    args = parser.parse_args()

    if args.save:
        matplotlib.use('Agg')

    args.h0, args.u0, source = base_flow(P, args.h0, args.u0, args.from_output)
    print(f"Base flow: {source}")
    Fr0 = args.u0 / np.sqrt(g*args.h0)
    c0 = np.sqrt(g*args.h0)
    wavelengths = np.linspace(args.wl_min, args.wl_max, args.n)
    D_bed = bed_diffusivity_from_bed_alpha(args.bed_alpha, args.dx)

    print(f"Base state: h0={args.h0:.5f} m, u0={args.u0:.4f} m/s, Fr0={Fr0:.4f}, "
          f"c0=sqrt(g h0)={c0*1000:.1f} mm/s")
    print(f"(settling velocity ws_g = {WS_G*1000:.3f} mm/s)")
    if args.bed_alpha:
        print(f"bed_alpha={args.bed_alpha} at dx={args.dx} m  ->  D_bed={D_bed:.4e} m^2/s")

    if args.sweep_kappa:
        mu_values = [0.01, 0.3, 0.6, 0.85, 0.95, 0.99, 0.995]
        results = [run_one(args.h0, args.u0, mu, args.cd_dep, wavelengths, sgn=args.sgn) for mu in mu_values]

        print(f"\n{'='*78}\nSummary across mu_over_sigma (kappa decreasing as "
              f"mu_over_sigma -> 1):\n{'='*78}")
        any_growth = False
        for mu, r in zip(mu_values, results):
            i_peak = np.argmax(r['bed_growth'])
            peak = r['bed_growth'][i_peak]
            any_growth = any_growth or peak > 0
            print(f"  mu_over_sigma={mu:6.3f}  kappa={r['kappa']:9.5f} 1/s  "
                  f"max bed-branch growth={peak:10.5f} 1/s at wl={wavelengths[i_peak]*100:6.2f} cm")
        if not any_growth:
            print("\n  -> NO mu_over_sigma in (0,1) (i.e. no kappa in this physically allowed "
                  "range) makes the bed/Exner branch unstable at this (h0,u0). The lag alone, "
                  "within the range the stochastic closure allows, is not sufficient here.")
        else:
            print("\n  -> At least one mu_over_sigma value DOES cross into instability - see above.")

        fig, axs = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True)
        cmap = plt.cm.viridis(np.linspace(0, 1, len(mu_values)))
        axs[0].axhline(0, color='gray', linewidth=0.7)
        for mu, r, c in zip(mu_values, results, cmap):
            axs[0].plot(wavelengths*100, r['bed_growth'], color=c, linewidth=1.8,
                        label=f'mu/sigma={mu} (kappa={r["kappa"]:.3f}/s)')
        axs[0].set_ylabel('bed/Exner branch growth rate [1/s]')
        axs[0].set_title(f'4-eq dispersion relation vs. kappa (lag): h0={args.h0:.4f}m, '
                          f'u0={args.u0:.3f}m/s, Fr0={Fr0:.3f}')
        axs[0].axvspan(11, 12, color='green', alpha=0.15, label='target 11-12cm')
        axs[0].legend(fontsize=7, ncol=2)

        for mu, r, c in zip(mu_values, results, cmap):
            axs[1].plot(wavelengths*100, r['bed_speed']*1000, color=c, linewidth=1.8)
        axs[1].axhline(0, color='gray', linewidth=0.7)
        axs[1].axvspan(11, 12, color='green', alpha=0.15)
        axs[1].set_xlabel('wavelength (cm)')
        axs[1].set_ylabel('migration speed (mm/s)\n(negative = upstream)')
        fig.tight_layout()

    elif args.sweep_bed_alpha:
        ba_values = [0.0, 0.003, 0.005, 0.0075, 0.008, 0.01]
        results = [run_one(args.h0, args.u0, args.mu_over_sigma, args.cd_dep, wavelengths,
                            D_bed=bed_diffusivity_from_bed_alpha(ba, args.dx), sgn=args.sgn)
                   for ba in ba_values]

        print(f"\n{'='*78}\nSummary across bed_alpha (dx={args.dx} m, "
              f"mu_over_sigma={args.mu_over_sigma}, cd_dep={args.cd_dep}):\n{'='*78}")
        for ba, r in zip(ba_values, results):
            i_peak = np.argmax(r['bed_growth'])
            print(f"  bed_alpha={ba:7.4f}  (D_bed={bed_diffusivity_from_bed_alpha(ba, args.dx):.3e} m^2/s)  "
                  f"peak growth={r['bed_growth'][i_peak]:10.5f} 1/s at wl={wavelengths[i_peak]*100:6.2f} cm, "
                  f"speed={r['bed_speed'][i_peak]*1000:7.3f} mm/s")

        fig, axs = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True)
        cmap = plt.cm.plasma(np.linspace(0, 0.85, len(ba_values)))
        axs[0].axhline(0, color='gray', linewidth=0.7)
        for ba, r, c in zip(ba_values, results, cmap):
            axs[0].plot(wavelengths*100, r['bed_growth'], color=c, linewidth=1.8,
                        label=f'bed_alpha={ba}')
        axs[0].set_ylabel('bed/Exner branch growth rate [1/s]')
        axs[0].set_title(f'4-eq dispersion relation vs. bed diffusion: h0={args.h0:.4f}m, '
                          f'u0={args.u0:.3f}m/s, mu/sigma={args.mu_over_sigma}, kappa fixed')
        axs[0].axvspan(11, 12, color='green', alpha=0.15, label='target 11-12cm')
        axs[0].legend(fontsize=8)

        for ba, r, c in zip(ba_values, results, cmap):
            axs[1].plot(wavelengths*100, r['bed_speed']*1000, color=c, linewidth=1.8)
        axs[1].axhline(0, color='gray', linewidth=0.7)
        axs[1].axvspan(11, 12, color='green', alpha=0.15)
        axs[1].set_xlabel('wavelength (cm)')
        axs[1].set_ylabel('migration speed (mm/s)\n(negative = upstream)')
        fig.tight_layout()

    else:
        r = run_one(args.h0, args.u0, args.mu_over_sigma, args.cd_dep, wavelengths, D_bed=D_bed,
                    sgn=args.sgn)

        fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        axs[0].axhline(0, color='gray', linewidth=0.7)
        for b in range(4):
            axs[0].plot(wavelengths*100, r['growth'][:, b], color='lightgray', linewidth=1)
        axs[0].plot(wavelengths*100, r['bed_growth'], color='crimson', linewidth=2,
                    label='bed/Exner-like branch')
        axs[0].plot(wavelengths*100, r['act_growth'], color='steelblue', linewidth=1.5,
                    linestyle='--', label='activity-like branch')
        axs[0].set_ylabel('growth rate Re(sigma)  [1/s]')
        axs[0].set_title(f'4-eq dispersion relation: h0={args.h0:.4f}m, u0={args.u0:.3f}m/s, '
                          f'Fr0={Fr0:.3f}, mu/sigma={args.mu_over_sigma}, cd_dep={args.cd_dep}, '
                          f'kappa={r["kappa"]:.4f}/s, bed_alpha={args.bed_alpha}')
        axs[0].axvspan(11, 12, color='green', alpha=0.15, label='target 11-12cm')
        axs[0].legend(fontsize=8)

        for b in range(4):
            axs[1].plot(wavelengths*100, r['speed'][:, b]*1000, color='lightgray', linewidth=1)
        axs[1].plot(wavelengths*100, r['bed_speed']*1000, color='crimson', linewidth=2)
        axs[1].axhline(0, color='gray', linewidth=0.7)
        axs[1].axvspan(11, 12, color='green', alpha=0.15)
        axs[1].set_xlabel('wavelength (cm)')
        axs[1].set_ylabel('migration speed (mm/s)\n(negative = upstream)')
        fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"Saved to {args.save}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
