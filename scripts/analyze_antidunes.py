#!/usr/bin/env python3
"""
analyze_antidunes.py - Hovmoeller (x-t) diagram of the band-pass-filtered
bed-thickness anomaly, with the dominant wavelength, RMS amplitude and
migration speed (positive = downstream, negative = upstream), and the
mean wavelength/amplitude of the antidunes: bedforms between upward zero
crossings with amplitude >= --min-amplitude (default 0, all).
Default wavelength band: 5-20 cm.

Usage:
    python scripts/analyze_antidunes.py OUTDIR [--xmin X] [--xmax X] [--band LO HI]
        [--min-amplitude A] [--save FILE.png]

Example:
    python scripts/analyze_antidunes.py _output --save antidunes.png
"""
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from plot_results import available_frames, read_frame, read_zbase


def bandpass(signal, dx, lo, hi):
    """Zero out every Fourier component outside wavelength [lo, hi] (m)."""
    n = len(signal)
    freq = np.fft.rfftfreq(n, d=dx)          # cycles/m
    spec = np.fft.rfft(signal)
    wavelength = np.full_like(freq, np.inf)
    wavelength[1:] = 1.0 / freq[1:]
    keep = (wavelength >= lo) & (wavelength <= hi)
    spec_filtered = np.where(keep, spec, 0.0)
    return np.fft.irfft(spec_filtered, n=n), freq, spec


def dominant_wavelength(signal, dx, lo, hi):
    n = len(signal)
    freq = np.fft.rfftfreq(n, d=dx)
    spec = np.abs(np.fft.rfft(signal))
    wavelength = np.full_like(freq, np.inf)
    wavelength[1:] = 1.0 / freq[1:]
    in_band = (wavelength >= lo) & (wavelength <= hi)
    if not np.any(in_band):
        return np.nan, 0.0
    idx = np.where(in_band)[0]
    peak = idx[np.argmax(spec[idx])]
    return wavelength[peak], spec[peak]


def bedform_statistics(x, thickness, lo, hi, min_amplitude=0.0):
    """Mean wavelength and amplitude of the antidunes in one profile.

    The linearly detrended thickness is band-passed to [lo, hi]; each
    bedform runs between two successive upward zero crossings (located by
    linear interpolation). Its wavelength is the distance between them,
    its amplitude half its crest-to-trough height. Only bedforms with
    amplitude >= min_amplitude count as antidunes.
    Returns (mean wavelength, mean amplitude, number of antidunes,
    RMS of the band-passed profile)."""
    dx = x[1] - x[0]
    anomaly = thickness - np.polyval(np.polyfit(x, thickness, 1), x)
    eta, *_ = bandpass(anomaly, dx, lo, hi)
    rms = float(np.sqrt(np.mean(eta**2)))
    up = np.where((eta[:-1] < 0.0) & (eta[1:] >= 0.0))[0]
    if len(up) < 2:
        return np.nan, np.nan, 0, rms
    x_cross = x[up] + dx * (-eta[up]) / (eta[up+1] - eta[up])
    wavelengths = np.diff(x_cross)
    amplitudes = np.array([0.5 * (eta[i0:i1+1].max() - eta[i0:i1+1].min())
                           for i0, i1 in zip(up[:-1], up[1:])])
    keep = amplitudes >= min_amplitude
    if not np.any(keep):
        return np.nan, np.nan, 0, rms
    return wavelengths[keep].mean(), float(amplitudes[keep].mean()), int(keep.sum()), rms


def migration_speed(sig_early, sig_late, dt, dx, lo, hi):
    """Speed = shift/dt, with the shift from the cross-correlation of two
    band-passed profiles. Positive = downstream, negative = upstream."""
    a, *_ = bandpass(sig_early, dx, lo, hi)
    b, *_ = bandpass(sig_late, dx, lo, hi)
    n = len(a)
    corr = np.fft.irfft(np.fft.rfft(b) * np.conj(np.fft.rfft(a)), n=n)
    # corr[k] = correlation at shift k*dx (circular); use the peak nearest 0
    # by looking only at shifts within +/- n/2.
    shifts = np.fft.fftfreq(n, d=1.0/n).astype(int)  # -n/2..n/2 in index units
    order = np.argsort(shifts)
    corr_ordered = corr[order]
    shifts_ordered = shifts[order]
    best = shifts_ordered[np.argmax(corr_ordered)]
    shift_m = best * dx
    return shift_m / dt


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('outdir')
    parser.add_argument('--xmin', type=float, default=None,
                         help='analysis window start (default: 25%% into the domain)')
    parser.add_argument('--xmax', type=float, default=None,
                         help='analysis window end (default: 90%% into the domain)')
    parser.add_argument('--band', type=float, nargs=2, default=[0.05, 0.20],
                         metavar=('LO', 'HI'),
                         help='wavelength band to search, meters (default 0.05 0.20)')
    parser.add_argument('--min-amplitude', type=float, default=0.0,
                         help='smallest amplitude counted as an antidune, meters (default 0)')
    parser.add_argument('--save', metavar='FILE.png')
    args = parser.parse_args()

    if args.save:
        matplotlib.use('Agg')

    frames = available_frames(args.outdir)
    z_base = read_zbase(args.outdir)
    lo, hi = args.band

    xs, ts, thick_rows = None, [], []
    for n in frames:
        x, h, hu, z, b, eta, t = read_frame(args.outdir, n)
        if xs is None:
            xs = x
            xmin = args.xmin if args.xmin is not None else x[0] + 0.25*(x[-1]-x[0])
            xmax = args.xmax if args.xmax is not None else x[0] + 0.90*(x[-1]-x[0])
            mask = (x >= xmin) & (x <= xmax)
            dx = x[1] - x[0]
        thickness = z - z_base
        thick_rows.append(thickness[mask])
        ts.append(t)
    ts = np.array(ts)
    xw = xs[mask]
    H = np.array(thick_rows)   # shape (n_times, n_x)

    # Band-pass each time row for the Hovmoeller diagram.
    H_filt = np.array([bandpass(row, dx, lo, hi)[0] for row in H])

    # Dominant wavelength averaged over the second half of the run (skip
    # the early transient).
    half = len(ts) // 2
    wavelengths = []
    for row in H[half:]:
        row_detrended = row - row.mean()
        wl, _ = dominant_wavelength(row_detrended, dx, lo, hi)
        if not np.isnan(wl):
            wavelengths.append(wl)
    mean_wavelength = np.mean(wavelengths) if wavelengths else float('nan')

    # Amplitude trend: RMS of the filtered signal at several times.
    rms = np.sqrt(np.mean(H_filt**2, axis=1))

    # Migration speed: compare the band-passed field at the start and end
    # of the second half of the run (avoids the initial transient).
    speed = migration_speed(H[half] - H[half].mean(), H[-1] - H[-1].mean(),
                             ts[-1]-ts[half], dx, lo, hi)

    print(f"Analysis window: x in [{xw[0]:.3f}, {xw[-1]:.3f}] m, "
          f"wavelength band [{lo:.3f}, {hi:.3f}] m")
    print(f"Dominant wavelength (2nd half of run, mean over frames): "
          f"{mean_wavelength*100:.2f} cm" if not np.isnan(mean_wavelength) else
          "Dominant wavelength: no clear peak in band")
    bedforms = [bedform_statistics(xw, row, lo, hi, args.min_amplitude) for row in H[half:]]
    if any(bf[2] > 0 for bf in bedforms):
        bf_wl = np.nanmean([bf[0] for bf in bedforms])
        bf_amp = np.nanmean([bf[1] for bf in bedforms])
        print(f"Mean antidune wavelength / amplitude (amplitude >= "
              f"{args.min_amplitude*1000:.2f} mm, 2nd half of run): "
              f"{bf_wl*100:.2f} cm / {bf_amp*1000:.3f} mm")
    else:
        print(f"No antidunes with amplitude >= {args.min_amplitude*1000:.2f} mm "
              f"in the 2nd half of the run")
    print(f"Band-passed RMS amplitude: start={rms[0]:.3e} m  "
          f"mid={rms[len(rms)//2]:.3e} m  end={rms[-1]:.3e} m")
    print(f"Estimated migration speed (2nd half of run): {speed*1000:.4f} mm/s "
          f"({'UPSTREAM' if speed < 0 else 'downstream'})")

    fig, axs = plt.subplots(2, 1, figsize=(9, 8))
    im = axs[0].imshow(H_filt, aspect='auto', origin='lower',
                        extent=[xw[0], xw[-1], ts[0], ts[-1]], cmap='RdBu_r')
    axs[0].set_xlabel('x (m)')
    axs[0].set_ylabel('t (s)')
    axs[0].set_title(f'Band-passed bed-thickness anomaly ({lo*100:.0f}-{hi*100:.0f} cm), '
                      f'speed={speed*1000:.3f} mm/s')
    fig.colorbar(im, ax=axs[0], label='thickness anomaly (m)')

    axs[1].plot(ts, rms)
    axs[1].set_xlabel('t (s)')
    axs[1].set_ylabel('band-passed RMS amplitude (m)')
    axs[1].set_title('Amplitude over time (growth check)')
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f'Saved to {args.save}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
