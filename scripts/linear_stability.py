#!/usr/bin/env python3
"""
linear_stability.py - linear growth rate and migration speed versus
wavelength for the 3-equation shallow-water/Exner system with the
equilibrium closure qs = qs_eq(h,u).

Linearises h_t+(hu)_x=0, (hu)_t+(hu^2+gh^2/2)_x=-gh(z_x+Sf) and
z_t+qs_x/(1-porosity)=0 around a uniform normal flow (h0,u0) (friction
slope = S0) for perturbations ~exp(i*k*x+sigma*t). The eigenvalues of
the 3x3 matrix L(k) give
  growth rate     = Re(sigma)      [1/s]
  migration speed = -Im(sigma)/k   [m/s]  (negative = upstream)

Unlike rp1.f90's cubic, this includes friction, the bed-slope source and
Qh = d(qs_eq)/dh.

All constants are read from the solver (model_params.py). The base flow
is the normal flow for Q_in unless --from-output (mean late-time flow of
a run) or --h0/--u0 are given.

Usage:
    python scripts/linear_stability.py [--h0 H --u0 U | --from-output DIR] [--save FILE.png]
"""
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from model_params import load, base_flow

P = load()
g, mann_n, d50, s_rel = P['g'], P['mann_n'], P['d50'], P['s_rel']
porosity, S0, tau_star_max, eps = P['porosity'], P['S0'], P['tau_star_max'], P['eps']
guo_c1, guo_p1, guo_c2, guo_p2, guo_p3 = (P[k] for k in ('guo_c1', 'guo_p1', 'guo_c2',
                                                          'guo_p2', 'guo_p3'))


def qs_equilibrium(h, u):
    """Guo/Recking equilibrium bedload law, as sediment_module.f90's
    qs_equilibrium."""
    cf = g * mann_n**2 / (abs(h)**(1.0/3.0) + eps)
    tau_star = min(cf * u**2 / (s_rel*g*d50 + eps), tau_star_max)
    Phi = (guo_c1*tau_star)**guo_p1 * (1.0 + (tau_star/guo_c2)**guo_p2)**(-guo_p3)
    return Phi * np.sqrt(s_rel*g*d50**3)


def Qu_Qh(h0, u0):
    eps_u = max(1.0e-4*max(abs(u0), 0.01), 1.0e-6)
    Qu = (qs_equilibrium(h0, u0+eps_u) - qs_equilibrium(h0, u0-eps_u)) / (2*eps_u)
    eps_h = max(1.0e-4*h0, 1.0e-8)
    Qh = (qs_equilibrium(h0+eps_h, u0) - qs_equilibrium(h0-eps_h, u0)) / (2*eps_h)
    return Qu, Qh


def L_matrix(k, h0, u0, Qu, Qh):
    return np.array([
        [-1j*k*u0,                        -1j*k*h0,              0.0],
        [-1j*k*g + g*4*S0/(3*h0),         -1j*k*u0 - g*2*S0/u0,  -1j*k*g],
        [-1j*k*Qh/(1-porosity),           -1j*k*Qu/(1-porosity), 0.0],
    ], dtype=complex)


def dispersion_curve(h0, u0, wavelengths):
    """All 3 eigenvalue branches, sorted by migration speed at each k."""
    Qu, Qh = Qu_Qh(h0, u0)
    n = len(wavelengths)
    growth = np.zeros((n, 3))
    speed = np.zeros((n, 3))
    for i, wl in enumerate(wavelengths):
        k = 2*np.pi / wl
        eigvals = np.linalg.eigvals(L_matrix(k, h0, u0, Qu, Qh))
        s = -eigvals.imag / k
        order = np.argsort(s)   # consistent ordering: slow/upstream -> fast/downstream
        growth[i] = eigvals.real[order]
        speed[i] = s[order]
    return growth, speed, Qu, Qh


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--h0', type=float, help='base flow depth (m)')
    parser.add_argument('--u0', type=float, help='base flow velocity (m/s)')
    parser.add_argument('--from-output', metavar='DIR',
                         help='use the mean flow of the 2nd half of this run as base flow')
    parser.add_argument('--wl-min', type=float, default=0.02)
    parser.add_argument('--wl-max', type=float, default=0.60)
    parser.add_argument('--n', type=int, default=400)
    parser.add_argument('--save', metavar='FILE.png')
    args = parser.parse_args()

    if args.save:
        matplotlib.use('Agg')

    args.h0, args.u0, source = base_flow(P, args.h0, args.u0, args.from_output)
    print(f"Base flow: {source}")
    Fr0 = args.u0 / np.sqrt(g*args.h0)
    c0 = np.sqrt(g*args.h0)
    wavelengths = np.linspace(args.wl_min, args.wl_max, args.n)
    growth, speed, Qu, Qh = dispersion_curve(args.h0, args.u0, wavelengths)

    print(f"Base state: h0={args.h0:.5f} m, u0={args.u0:.4f} m/s, Fr0={Fr0:.4f}, "
          f"c0=sqrt(g h0)={c0*1000:.1f} mm/s")
    print(f"Qu = d(qs_eq)/du = {Qu:.4e},  Qh = d(qs_eq)/dh = {Qh:.4e}")
    print(f"(for reference, pure gravity-wave speeds would be u0-c0={1000*(args.u0-c0):.1f} "
          f"and u0+c0={1000*(args.u0+c0):.1f} mm/s - the two branches here that sit near "
          f"those values are ordinary water waves, not bed-related; the third, slow branch "
          f"is the bed/Exner-like one relevant to antidunes)")

    # Bed branch: the root whose speed is furthest from both u0 - c0 and u0 + c0
    gravity_speeds = np.array([args.u0 - c0, args.u0 + c0])
    dist_to_gravity = np.min(np.abs(speed[:, :, None] - gravity_speeds[None, None, :]), axis=2)
    bed_branch = np.argmax(dist_to_gravity, axis=1)
    bed_growth = growth[np.arange(len(wavelengths)), bed_branch]
    bed_speed = speed[np.arange(len(wavelengths)), bed_branch]

    if np.max(bed_growth) > 0:
        i_peak = np.argmax(bed_growth)
        print(f"\nBed/Exner branch: fastest-growing wavelength "
              f"{wavelengths[i_peak]*100:.2f} cm (growth rate {bed_growth[i_peak]:.4f} 1/s, "
              f"migration speed {bed_speed[i_peak]*1000:.3f} mm/s "
              f"{'UPSTREAM' if bed_speed[i_peak]<0 else 'downstream'})")
    else:
        i_peak = np.argmax(bed_growth)
        print(f"\nBed/Exner branch: NO growing wavelength in this band - least-decaying "
              f"is {bed_growth[i_peak]:.4e} 1/s at {wavelengths[i_peak]*100:.2f} cm "
              f"(migration speed there: {bed_speed[i_peak]*1000:.3f} mm/s)")

    fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axs[0].axhline(0, color='gray', linewidth=0.7)
    for b in range(3):
        axs[0].plot(wavelengths*100, growth[:, b], color='lightgray', linewidth=1,
                    label='other branches' if b == 0 else None)
    axs[0].plot(wavelengths*100, bed_growth, color='crimson', linewidth=2,
                label='bed/Exner-like branch')
    axs[0].set_ylabel('growth rate Re(sigma)  [1/s]')
    axs[0].set_title(f'Linear dispersion relation: h0={args.h0:.4f}m, '
                      f'u0={args.u0:.3f}m/s, Fr0={Fr0:.3f}')
    axs[0].axvspan(8, 15, color='green', alpha=0.15, label='observed antidune band (8-15cm)')
    axs[0].legend(fontsize=8)

    for b in range(3):
        axs[1].plot(wavelengths*100, speed[:, b]*1000, color='lightgray', linewidth=1)
    axs[1].plot(wavelengths*100, bed_speed*1000, color='crimson', linewidth=2)
    axs[1].axhline(0, color='gray', linewidth=0.7)
    axs[1].axvspan(8, 15, color='green', alpha=0.15)
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
