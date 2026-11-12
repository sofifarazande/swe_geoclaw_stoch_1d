#!/usr/bin/env python3
"""
plot_results.py - plot the free surface eta = h+z, the bed z, velocity u
and particle activity b from fort.qNNNN / fort.tNNNN files (numpy and
matplotlib only). Works for the main case and both tests.

Usage:
    python scripts/plot_results.py OUTDIR [--frame N] [--animate] [--save FILE]

Examples:
    python scripts/plot_results.py _output                    # last frame, show on screen
    python scripts/plot_results.py _output --frame 50          # a specific frame
    python scripts/plot_results.py _output --animate                      # animate, show on screen
    python scripts/plot_results.py _output --animate --save run.mp4       # animate, save to file
    python scripts/plot_results.py tests/dam_break/_output --frame -1 --save last_frame.png
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib

# Non-GUI backend when saving; must be chosen before importing pyplot.
if '--save' in sys.argv:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.animation as animation


def available_frames(outdir):
    """Sorted list of frame numbers found in outdir (from fort.qNNNN files)."""
    files = glob.glob(os.path.join(outdir, 'fort.q[0-9][0-9][0-9][0-9]'))
    frames = sorted(int(re.search(r'fort\.q(\d{4})', f).group(1)) for f in files)
    if not frames:
        raise FileNotFoundError(f'No fort.qNNNN files found in {outdir!r}')
    return frames


def read_frame(outdir, n):
    """Read one output frame: returns x, h, hu, z, b, eta, t (all 1D arrays except t)."""
    qfile = os.path.join(outdir, f'fort.q{n:04d}')
    tfile = os.path.join(outdir, f'fort.t{n:04d}')

    with open(qfile) as f:
        f.readline()                      # grid_number
        f.readline()                      # AMR_level
        mx = int(f.readline().split()[0])
        xlow = float(f.readline().split()[0])
        dx = float(f.readline().split()[0])

    data = np.loadtxt(qfile, skiprows=6, max_rows=mx)
    h, hu, z, b, eta = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    x = xlow + (np.arange(mx) + 0.5) * dx

    with open(tfile) as f:
        t = float(f.readline().split()[0])

    return x, h, hu, z, b, eta, t


def read_zbase(outdir, n=0):
    """Read the rigid floor z_base = aux(1,:) from fort.aNNNN (written
    once, at frame 0)."""
    afile = os.path.join(outdir, f'fort.a{n:04d}')
    with open(afile) as f:
        f.readline(); f.readline()        # grid_number, AMR_level
        mx = int(f.readline().split()[0])
    data = np.loadtxt(afile, skiprows=6, max_rows=mx)
    return data[:, 0]   # aux(1,:) = z_base


def velocity(h, hu, dry_tol=1.e-6):
    u = np.zeros_like(h)
    wet = h > dry_tol
    u[wet] = hu[wet] / h[wet]
    return u


def make_axes():
    fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    axs[0].set_ylabel('elevation (m)')
    axs[1].set_ylabel('u (m/s)')
    axs[2].set_ylabel('b (particle activity)')
    axs[2].set_xlabel('x (m)')
    return fig, axs


def plot_frame(outdir, n, axs=None, fig=None):
    x, h, hu, z, b, eta, t = read_frame(outdir, n)
    u = velocity(h, hu)

    if axs is None:
        fig, axs = make_axes()

    for ax in axs:
        ax.cla()

    axs[0].plot(x, eta, 'b-', label='eta = h+z')
    axs[0].plot(x, z, 'k-', label='z (bed)')
    axs[0].legend(loc='upper right')
    axs[0].set_ylabel('elevation (m)')

    axs[1].plot(x, u, 'r-')
    axs[1].set_ylabel('u (m/s)')

    axs[2].plot(x, b, 'g-')
    axs[2].set_ylabel('b')
    axs[2].set_xlabel('x (m)')

    fig.suptitle(f'frame {n},  t = {t:.4f} s')
    return fig, axs


def animate(outdir, save=None):
    frames = available_frames(outdir)
    fig, axs = make_axes()

    def update(n):
        plot_frame(outdir, n, axs=axs, fig=fig)
        print(f'  encoding frame {n} ({frames.index(n)+1}/{len(frames)})', end='\r')

    ani = animation.FuncAnimation(fig, update, frames=frames, interval=120)

    if save:
        print(f'Encoding {len(frames)} frames to {save} - this can take a while '
              'with no other output in the meantime, that is normal...')
        try:
            ani.save(save)
        except Exception as e:
            print(f'\nFailed to save {save}: {e}')
            if save.endswith('.mp4'):
                print('(.mp4 needs the ffmpeg command available on PATH - '
                      'try --save with a .gif file instead, which only needs '
                      'matplotlib\'s own built-in Pillow writer.)')
            raise
        print(f'\nSaved animation to {save}')
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('outdir', help='directory with fort.qNNNN/fort.tNNNN files')
    parser.add_argument('--frame', type=int, default=-1,
                         help='frame number to plot (default: last available); '
                              'negative indices count from the end, like a Python list')
    parser.add_argument('--animate', action='store_true',
                         help='animate over all available frames instead of one')
    parser.add_argument('--save', metavar='FILE',
                         help='save to FILE (.png for a single frame, .mp4/.gif for --animate) '
                              'instead of showing on screen')
    args = parser.parse_args()

    if args.animate:
        animate(args.outdir, save=args.save)
    else:
        frames = available_frames(args.outdir)
        n = frames[args.frame] if args.frame < 0 else args.frame
        plot_frame(args.outdir, n)
        if args.save:
            plt.savefig(args.save, dpi=150)
            print(f'Saved to {args.save}')
        else:
            plt.show()


if __name__ == '__main__':
    main()
