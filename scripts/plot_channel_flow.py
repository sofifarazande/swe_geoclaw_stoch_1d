#!/usr/bin/env python3
"""
plot_channel_flow.py - side view of the channel: rigid floor z_base
(gray), mobile sediment layer z (brown) and water (blue, wet cells only).
Uses the frame readers in plot_results.py.

Usage:
    python scripts/plot_channel_flow.py OUTDIR [--frame N] [--animate] [--save FILE]

Examples:
    python scripts/plot_channel_flow.py _output --animate                  # show on screen
    python scripts/plot_channel_flow.py _output --animate --save flow.mp4  # save a video
    python scripts/plot_channel_flow.py _output --frame -1 --save last.png # last frame only
"""
import argparse
import sys

import numpy as np
import matplotlib

if '--save' in sys.argv:
    matplotlib.use('Agg')   # non-GUI backend, before importing pyplot

import matplotlib.pyplot as plt
import matplotlib.animation as animation

from plot_results import available_frames, read_frame, read_zbase


def channel_bounds(outdir, frames, z_base):
    """Axis limits over all frames, so the y-axis stays fixed during an
    animation."""
    z_min, z_max, eta_max = np.inf, -np.inf, -np.inf
    for n in frames:
        x, h, hu, z, b, eta, t = read_frame(outdir, n)
        z_min = min(z_min, z.min())
        z_max = max(z_max, z.max())
        eta_max = max(eta_max, eta.max())
    z_base_min = z_base.min()
    return min(z_min, z_base_min), z_max, eta_max, z_base_min


def plot_channel(outdir, n, ax, bounds, z_base, dry_tol=1.e-6):
    x, h, hu, z, b, eta, t = read_frame(outdir, n)
    wet = h > dry_tol
    z_min, z_max, eta_max, z_base_min = bounds

    ax.cla()

    # Margin relative to the total plot height, so a flat bed still shows the ground band
    ground_bottom = z_base_min - 0.15 * (eta_max - z_min + 1.e-3)

    # Rigid floor z_base
    ax.fill_between(x, ground_bottom, z_base, color='dimgray', zorder=0)
    ax.plot(x, z_base, color='k', linewidth=1.2, zorder=1)

    # Mobile sediment layer between z_base and z
    ax.fill_between(x, z_base, z, color='saddlebrown', zorder=2)
    ax.plot(x, z, color='k', linewidth=1, zorder=3)                 # sediment-surface line

    eta_wet = np.where(wet, eta, np.nan)
    ax.fill_between(x, z, eta_wet, color='deepskyblue', alpha=0.85, zorder=4)
    ax.plot(x, eta_wet, color='navy', linewidth=1, zorder=5)        # free-surface line

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(ground_bottom, eta_max + 0.08 * (eta_max - z_min + 1.e-3))
    ax.set_xlabel('x (m)')
    ax.set_ylabel('elevation (m)')
    ax.set_title(f't = {t:.4f} s')


def animate(outdir, save=None):
    frames = available_frames(outdir)
    z_base = read_zbase(outdir)
    bounds = channel_bounds(outdir, frames, z_base)
    fig, ax = plt.subplots(figsize=(10, 3.5))

    def update(n):
        plot_channel(outdir, n, ax, bounds, z_base)

    ani = animation.FuncAnimation(fig, update, frames=frames, interval=120)

    if save:
        print(f'Encoding {len(frames)} frames to {save} - this can take a while '
              'with no other output in the meantime, that is normal...')
        try:
            ani.save(save)
        except Exception as e:
            print(f'Failed to save {save}: {e}')
            if save.endswith('.mp4'):
                print("(.mp4 needs ffmpeg on PATH - try a .gif file instead, "
                      "which only needs matplotlib's bundled Pillow writer.)")
            raise
        print(f'Saved animation to {save}')
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
        z_base = read_zbase(args.outdir)
        bounds = channel_bounds(args.outdir, frames, z_base)
        fig, ax = plt.subplots(figsize=(10, 3.5))
        plot_channel(args.outdir, n, ax, bounds, z_base)
        if args.save:
            plt.savefig(args.save, dpi=150)
            print(f'Saved to {args.save}')
        else:
            plt.show()


if __name__ == '__main__':
    main()
