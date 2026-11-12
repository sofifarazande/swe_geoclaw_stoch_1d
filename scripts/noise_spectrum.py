#!/usr/bin/env python3
"""
noise_spectrum.py - stationary, noise-driven power spectrum of the
linearised 4-equation system of linear_stability_4eq.py, with noise
entering b only, as in src1.f90.

For dv/dt = L(k) v + noise, v = (h',u',z',b'), the stationary covariance
C(k) solves the Lyapunov equation

    L(k) C(k) + C(k) L(k)^H = -Q(k),

with Q nonzero only in the (b,b) entry. C_zz(k) is the predicted bed
spectrum; its peak is the predicted dominant wavelength.

The block-averaged noise (one increment per block of length ell_block,
independent between blocks) has spectral density

    Q_bb(k) = amp^2 * ell_block * sinc^2(k*ell_block/2),

with amp^2 = 2*mu*b0*v_p_grain/(ell_block*B_channel), as in src1.f90.
Blocks are treated as randomly placed, and the Euler-Maruyama substeps
as the continuous-time limit.

All constants are read from the solver (model_params.py); run-time
values default to setprob.data/claw.data, and the base flow is chosen as
in linear_stability_4eq.py.

Usage:
    python scripts/noise_spectrum.py [--h0 H --u0 U | --from-output DIR]
        [--mu-over-sigma M] [--cd-dep C] [--bed-alpha A] [--n-diam-reaction N]
        [--save FILE.png]
    python scripts/noise_spectrum.py --compare-blocks [--measured W1 W2 W3] [--save FILE.png]
        (predicted S_z(k) for n_diam_reaction = 7.5/15/30, optionally with
        the wavelengths [cm] measured in the corresponding simulations)
"""
import argparse

import numpy as np
from scipy.linalg import solve_continuous_lyapunov
import matplotlib
import matplotlib.pyplot as plt

from linear_stability_4eq import (
    P, d50, kappa_of, base_state, derivatives, L_matrix,
    bed_diffusivity_from_bed_alpha,
)
from model_params import base_flow

B_channel = P['B_channel']
v_p_grain = P['v_p_grain']


def block_length(n_diam_reaction, dx):
    """ell_block = block_ncells*dx, as in init_stochastic_transport."""
    ell_reaction = n_diam_reaction * d50
    block_ncells = max(1, round(ell_reaction / dx))
    return block_ncells * dx


def noise_amp_squared(mu_g, b0, ell_block):
    """amp^2 = 2*mu*b0*v_p_grain/(ell_block*B_channel), as in src1.f90
    with the base-state b0."""
    k_factor = v_p_grain / (ell_block * B_channel)
    return 2.0 * mu_g * b0 * k_factor


def Q_bb(k, amp2, ell_block):
    """Spectral density of the block-averaged noise (module docstring)."""
    x = k * ell_block / 2.0
    sinc = np.sin(x)/x if abs(x) > 1e-12 else 1.0
    return amp2 * ell_block * sinc**2


def stationary_spectrum(h0, u0, mu_over_sigma, cd_dep, wavelengths, dx, n_diam_reaction,
                         D_bed=0.0):
    kappa, sigma_g, mu_g = kappa_of(mu_over_sigma, cd_dep)
    u_p0, lambda0, b0, Du0 = base_state(h0, u0, kappa)
    Up_h, Up_u, Up_S, Lam_h, Lam_u, Lam_S = derivatives(h0, u0, kappa)
    ell_block = block_length(n_diam_reaction, dx)
    amp2 = noise_amp_squared(mu_g, b0, ell_block)

    S_h = np.zeros(len(wavelengths))
    S_u = np.zeros(len(wavelengths))
    S_z = np.zeros(len(wavelengths))
    S_b = np.zeros(len(wavelengths))
    for i, wl in enumerate(wavelengths):
        k = 2*np.pi/wl
        L = L_matrix(k, h0, u0, kappa, u_p0, b0, Up_h, Up_u, Up_S, Lam_h, Lam_u, Lam_S, D_bed, Du0)
        Q = np.zeros((4, 4), dtype=complex)
        Q[3, 3] = Q_bb(k, amp2, ell_block)
        C = solve_continuous_lyapunov(L, -Q)
        S_h[i] = C[0, 0].real
        S_u[i] = C[1, 1].real
        S_z[i] = C[2, 2].real
        S_b[i] = C[3, 3].real
    return dict(S_h=S_h, S_u=S_u, S_z=S_z, S_b=S_b, kappa=kappa, b0=b0, ell_block=ell_block,
                amp2=amp2)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--h0', type=float, help='base flow depth (m)')
    parser.add_argument('--u0', type=float, help='base flow velocity (m/s)')
    parser.add_argument('--from-output', metavar='DIR',
                         help='use the mean flow of the 2nd half of this run as base flow')
    parser.add_argument('--mu-over-sigma', type=float, default=P['mu_over_sigma'])
    parser.add_argument('--cd-dep', type=float, default=P['cd_dep'])
    parser.add_argument('--bed-alpha', type=float, default=P['bed_alpha'])
    parser.add_argument('--dx', type=float, default=P['dx'])
    parser.add_argument('--n-diam-reaction', type=float, default=P['n_diam_reaction'])
    parser.add_argument('--wl-min', type=float, default=0.03)
    parser.add_argument('--wl-max', type=float, default=0.40)
    parser.add_argument('--n', type=int, default=400)
    parser.add_argument('--compare-blocks', action='store_true',
                         help='overlay predicted S_z(k) for n_diam_reaction=7.5/15/30')
    parser.add_argument('--measured', type=float, nargs=3, metavar=('W7.5', 'W15', 'W30'),
                         help='with --compare-blocks: measured wavelengths [cm] of runs '
                              'with n_diam_reaction=7.5/15/30')
    parser.add_argument('--save', metavar='FILE.png')
    args = parser.parse_args()

    if args.save:
        matplotlib.use('Agg')

    args.h0, args.u0, source = base_flow(P, args.h0, args.u0, args.from_output)
    print(f"Base flow: {source}, h0={args.h0:.5f} m, u0={args.u0:.4f} m/s")
    wavelengths = np.linspace(args.wl_min, args.wl_max, args.n)
    D_bed = bed_diffusivity_from_bed_alpha(args.bed_alpha, args.dx)

    if args.compare_blocks:
        cases = list(zip([7.5, 15.0, 30.0], args.measured or [None]*3))
        fig, ax = plt.subplots(figsize=(9, 6))
        cmap = plt.cm.plasma(np.linspace(0, 0.8, len(cases)))
        print(f"{'='*78}\nPredicted vs. measured dominant wavelength, "
              f"bed_alpha={args.bed_alpha}, mu/sigma={args.mu_over_sigma}, "
              f"cd_dep={args.cd_dep}\n{'='*78}")
        for (ndr, measured_wl), c in zip(cases, cmap):
            r = stationary_spectrum(args.h0, args.u0, args.mu_over_sigma, args.cd_dep,
                                     wavelengths, args.dx, ndr, D_bed)
            i_peak = np.argmax(r['S_z'])
            predicted_wl = wavelengths[i_peak]*100
            ell_block = r['ell_block']
            line = (f"  n_diam_reaction={ndr:5.1f} (block={ell_block*100:.2f}cm): "
                    f"predicted peak={predicted_wl:6.2f}cm")
            if measured_wl:
                line += f"   measured={measured_wl:6.2f}cm   ratio={predicted_wl/measured_wl:.2f}"
            print(line)
            ax.plot(wavelengths*100, r['S_z']/r['S_z'].max(), color=c, linewidth=2,
                    label=f'n_diam_reaction={ndr} (block={ell_block*100:.1f}cm)')
            if measured_wl:
                ax.axvline(measured_wl, color=c, linestyle='--', alpha=0.6)
        ax.set_xlabel('wavelength (cm)')
        ax.set_ylabel('predicted S_z(k), normalized to peak=1')
        ax.set_title('Noise-driven stationary spectrum vs. block length'
                     + ('\n(dashed lines = measured wavelengths)' if args.measured else ''))
        ax.legend(fontsize=9)
        fig.tight_layout()

    else:
        r = stationary_spectrum(args.h0, args.u0, args.mu_over_sigma, args.cd_dep,
                                 wavelengths, args.dx, args.n_diam_reaction, D_bed)
        i_peak = np.argmax(r['S_z'])
        print(f"kappa={r['kappa']:.4f} 1/s, b0={r['b0']:.4e} m, "
              f"ell_block={r['ell_block']*100:.2f}cm, amp^2={r['amp2']:.4e}")
        print(f"Predicted dominant wavelength (S_z peak): {wavelengths[i_peak]*100:.2f} cm")

        fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        axs[0].plot(wavelengths*100, r['S_z'], color='crimson', linewidth=2, label='S_z (bed)')
        axs[0].set_ylabel('predicted stationary power spectrum')
        axs[0].axvspan(11, 12, color='green', alpha=0.15, label='target 11-12cm')
        axs[0].legend(fontsize=8)
        axs[0].set_title(f'Noise-driven stationary spectrum: h0={args.h0:.4f}m, '
                          f'u0={args.u0:.3f}m/s, mu/sigma={args.mu_over_sigma}, '
                          f'bed_alpha={args.bed_alpha}, n_diam_reaction={args.n_diam_reaction}')

        axs[1].plot(wavelengths*100, r['S_h'], label='S_h')
        axs[1].plot(wavelengths*100, r['S_u'], label='S_u')
        axs[1].plot(wavelengths*100, r['S_b'], label='S_b')
        axs[1].set_yscale('log')
        axs[1].set_xlabel('wavelength (cm)')
        axs[1].set_ylabel('other variables (log scale)')
        axs[1].legend(fontsize=8)
        fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"Saved to {args.save}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
