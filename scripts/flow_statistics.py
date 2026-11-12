#!/usr/bin/env python3
"""
flow_statistics.py - time series of domain statistics: mean depth,
mean/max velocity, bulk Froude number, wetted fraction, sediment
transport rate qs, sediment-layer thickness z - z_base, bed elevation z,
and the mean antidune wavelength and amplitude.

- mean_h, mean_u, max_u and Froude use wet cells only (h > dry_tol);
  Fr = |mean_u| / sqrt(g*mean_h).
- Thickness, z and qs are averaged over the whole domain; qs_outlet is
  qs in the last cell.
- qs = sign(u_p*b, u) is recomputed from h, u and b with the closure of
  sediment_module.f90, whose constants are read from the Fortran source
  (model_params.py).
- Antidunes: in the window [xmin, xmax], the thickness anomaly is
  band-passed to the wavelength band and split into bedforms at its
  upward zero crossings. Bedforms with amplitude (half the
  crest-to-trough height) >= --min-amplitude (default 0, all) are antidunes;
  antidune_wavelength and antidune_amplitude are their means (NaN if
  none), antidune_count their number (see
  analyze_antidunes.bedform_statistics). bed_rms is the RMS of the
  band-passed anomaly, with no threshold.

Usage:
    python scripts/flow_statistics.py OUTDIR [--save FILE.csv] [--dry-tol TOL]
        [--xmin X] [--xmax X] [--band LO HI] [--min-amplitude A]

Examples:
    python scripts/flow_statistics.py _output --save stats.csv
    python scripts/flow_statistics.py tests/dam_break/_output --save dam_break_stats.csv
"""
import argparse

import numpy as np

from plot_results import available_frames, read_frame, read_zbase
from analyze_antidunes import bedform_statistics
from model_params import load, settling_velocity

P = load()
g = P['g']
WS_G = settling_velocity(P)


def bed_slope(z, dx):
    """slope = -dz/dx: centred inside, one-sided at the ends."""
    slope = np.empty_like(z)
    slope[1:-1] = -(z[2:] - z[:-2]) / (2.0*dx)
    slope[0] = -(z[1] - z[0]) / dx
    slope[-1] = -(z[-1] - z[-2]) / dx
    return slope


def transport_rate(h, u, b, z, dx):
    """qs = sign(u_p*b, u), as qs_stochastic and particle_velocity in
    sediment_module.f90."""
    s_rel, d50, eps = P['s_rel'], P['d50'], P['eps']
    alpha_fric, min_us = P['alpha_fric'], P['min_us']
    slope = bed_slope(z, dx)

    cf = g * P['mann_n']**2 / (np.abs(h)**(1.0/3.0) + eps)
    tau_star = np.minimum(cf * u**2 / (s_rel * g * d50 + eps), P['tau_star_max'])
    theta = np.arctan(slope)
    tau_star_c = np.maximum(P['tau_star_0'] * np.sin(alpha_fric - theta) / np.sin(alpha_fric), 0.0)

    u_p = P['bagnold_b_coeff']*WS_G + P['bagnold_a']*P['sqrt_sgd'] * \
        (np.sqrt(np.maximum(tau_star, 0.0)) - np.sqrt(tau_star_c))
    u_p = np.maximum(u_p, min_us)
    u_p = np.where(h <= P['h_sed_min'], min_us, u_p)   # particle_velocity's low-h branch

    # Fortran sign(A,0) = +|A|, unlike np.sign(0) = 0
    qs = np.where(u >= 0.0, 1.0, -1.0) * u_p * b
    qs = np.where((h > P['h_wet']) & (b > 0.0), qs, 0.0)   # qs_stochastic's dry/no-activity gate
    return qs


def frame_statistics(outdir, n, z_base, window, band, min_amplitude, dry_tol=1.e-6):
    x, h, hu, z, b, eta, t = read_frame(outdir, n)
    wet = h > dry_tol
    dx = x[1] - x[0]

    n_wet = np.count_nonzero(wet)
    wetted_frac = n_wet / len(h)
    mean_z, min_z, max_z = z.mean(), z.min(), z.max()
    thickness = z - z_base
    mean_th, min_th, max_th = thickness.mean(), thickness.min(), thickness.max()

    # qs over the whole domain (zero in dry or inactive cells)
    u_all = np.zeros_like(h)
    u_all[wet] = hu[wet] / h[wet]
    qs = transport_rate(h, u_all, b, z, dx)
    mean_qs, max_qs = qs.mean(), np.max(np.abs(qs))
    # transport rate in the last cell (outlet)
    qs_outlet = qs[-1]

    in_window = (x >= window[0]) & (x <= window[1])
    ad_wl, ad_amp, ad_n, bed_rms = bedform_statistics(x[in_window], thickness[in_window],
                                                       *band, min_amplitude)

    common = dict(frame=n, t=t, wetted_frac=wetted_frac,
                  mean_th=mean_th, min_th=min_th, max_th=max_th,
                  mean_z=mean_z, min_z=min_z, max_z=max_z,
                  mean_qs=mean_qs, max_qs=max_qs, qs_outlet=qs_outlet,
                  ad_wl=ad_wl, ad_amp=ad_amp, ad_n=ad_n, bed_rms=bed_rms)

    if n_wet == 0:
        return dict(common, mean_h=0.0, mean_u=0.0, max_u=0.0, froude=np.nan)

    mean_h = h[wet].mean()
    u_wet = hu[wet] / h[wet]
    mean_u = u_wet.mean()
    max_u = np.max(np.abs(u_wet))
    froude = abs(mean_u) / np.sqrt(g * mean_h)

    return dict(common, mean_h=mean_h, mean_u=mean_u, max_u=max_u, froude=froude)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('outdir', help='directory with fort.qNNNN/fort.tNNNN files')
    parser.add_argument('--save', metavar='FILE.csv',
                         help='also write the table to a CSV file')
    parser.add_argument('--dry-tol', type=float, default=1.e-6,
                         help='depth below which a cell counts as dry (default 1e-6 m)')
    parser.add_argument('--xmin', type=float, default=None,
                         help='antidune window start (default: 25%% into the domain)')
    parser.add_argument('--xmax', type=float, default=None,
                         help='antidune window end (default: 90%% into the domain)')
    parser.add_argument('--band', type=float, nargs=2, default=[0.05, 0.20],
                         metavar=('LO', 'HI'),
                         help='antidune wavelength band, meters (default 0.05 0.20)')
    parser.add_argument('--min-amplitude', type=float, default=0.0,
                         help='smallest amplitude counted as an antidune, meters (default 0)')
    args = parser.parse_args()

    frames = available_frames(args.outdir)
    z_base = read_zbase(args.outdir)
    x = read_frame(args.outdir, frames[0])[0]
    window = (args.xmin if args.xmin is not None else x[0] + 0.25*(x[-1]-x[0]),
              args.xmax if args.xmax is not None else x[0] + 0.90*(x[-1]-x[0]))
    rows = [frame_statistics(args.outdir, n, z_base, window, args.band, args.min_amplitude,
                             args.dry_tol)
            for n in frames]

    print(f"Antidune window x in [{window[0]:.3f}, {window[1]:.3f}] m, "
          f"band [{args.band[0]:.3f}, {args.band[1]:.3f}] m, "
          f"amplitude >= {args.min_amplitude*1000:.2f} mm")
    header = (f"{'t':>10} {'mean_h':>10} {'mean_u':>10} {'max_u':>10} {'Froude':>10} "
              f"{'wetted%':>9} {'mean_qs':>11} {'max_qs':>11} {'qs_outlet':>11} "
              f"{'mean_th':>10} {'min_th':>10} {'max_th':>10} "
              f"{'mean_z':>10} {'min_z':>10} {'max_z':>10} "
              f"{'ad_wl_cm':>9} {'ad_amp_mm':>10} {'ad_n':>5} {'rms_mm':>8}")
    print(header)
    print('-' * len(header))
    for r in rows:
        print(f"{r['t']:10.4f} {r['mean_h']:10.5f} {r['mean_u']:10.5f} "
              f"{r['max_u']:10.5f} {r['froude']:10.4f} {100*r['wetted_frac']:9.1f} "
              f"{r['mean_qs']:11.3e} {r['max_qs']:11.3e} {r['qs_outlet']:11.3e} "
              f"{r['mean_th']:10.5f} {r['min_th']:10.5f} {r['max_th']:10.5f} "
              f"{r['mean_z']:10.5f} {r['min_z']:10.5f} {r['max_z']:10.5f} "
              f"{100*r['ad_wl']:9.2f} {1000*r['ad_amp']:10.4f} {r['ad_n']:5d} "
              f"{1000*r['bed_rms']:8.4f}")

    if args.save:
        with open(args.save, 'w') as f:
            f.write('frame,t,mean_h,mean_u,max_u,froude,wetted_frac,mean_qs,max_qs,qs_outlet,'
                    'mean_thickness,min_thickness,max_thickness,mean_z,min_z,max_z,'
                    'antidune_wavelength,antidune_amplitude,antidune_count,bed_rms\n')
            for r in rows:
                f.write(f"{r['frame']},{r['t']:.6f},{r['mean_h']:.6f},{r['mean_u']:.6f},"
                        f"{r['max_u']:.6f},{r['froude']:.6f},{r['wetted_frac']:.6f},"
                        f"{r['mean_qs']:.6e},{r['max_qs']:.6e},{r['qs_outlet']:.6e},"
                        f"{r['mean_th']:.6f},{r['min_th']:.6f},{r['max_th']:.6f},"
                        f"{r['mean_z']:.6f},{r['min_z']:.6f},{r['max_z']:.6f},"
                        f"{r['ad_wl']:.6f},{r['ad_amp']:.6e},{r['ad_n']},{r['bed_rms']:.6e}\n")
        print(f'\nSaved to {args.save}')


if __name__ == '__main__':
    main()
