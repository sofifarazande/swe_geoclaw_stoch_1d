"""
Compare the dam-break result at the final output time with Ritter's
analytical solution (flat bed, frictionless, dry downstream) and save an
overlay plot.

Ritter's solution, with xi = (x-x0)/t and c0 = sqrt(g*h0):
    xi <= -c0        : h = h0,               u = 0
    -c0 < xi < 2*c0  : h = (2c0-xi)^2/(9g),  u = (2/3)(xi+c0)
    xi >= 2*c0       : h = 0,                u = 0
x0 and h0 must match qinit.f90.
"""
import glob
import numpy as np

g = 9.81
x0 = 3.0
h0 = 0.5


def ritter(x, t):
    c0 = np.sqrt(g * h0)
    xi = (x - x0) / t
    h = np.where(xi <= -c0, h0, np.where(xi >= 2*c0, 0.0, (2*c0 - xi)**2 / (9.0*g)))
    u = np.where(xi <= -c0, 0.0, np.where(xi >= 2*c0, 0.0, (2.0/3.0)*(xi + c0)))
    return h, u


def main(outdir='_output'):
    qfiles = sorted(glob.glob(f'{outdir}/fort.q[0-9][0-9][0-9][0-9]'))
    if not qfiles:
        print(f'No fort.q files found in {outdir} - run the case first (make .output).')
        return
    fname = qfiles[-1]
    tname = fname.replace('fort.q', 'fort.t')
    with open(tname) as f:
        t = float(f.readline().split()[0])

    data = np.loadtxt(fname, skiprows=5)
    h_num, hu_num = data[:,0], data[:,1]
    with open(fname) as f:
        f.readline(); f.readline()
        mx = int(f.readline().split()[0])
        xlow = float(f.readline().split()[0])
        dx = float(f.readline().split()[0])
    x = xlow + (np.arange(mx) + 0.5) * dx

    u_num = np.zeros_like(h_num)
    wet = h_num > 1.e-6
    u_num[wet] = hu_num[wet] / h_num[wet]

    h_ex, u_ex = ritter(x, t)

    err_h = np.abs(h_num - h_ex)
    err_u = np.abs(u_num - u_ex)

    print(f"Comparing against Ritter's solution at t={t:.4f}s "
          f"(x0={x0}, h0={h0})")
    print(f"  L1 error in h:   {np.mean(err_h):.4e}   Linf error in h:   {np.max(err_h):.4e}")
    print(f"  L1 error in u:   {np.mean(err_u):.4e}   Linf error in u:   {np.max(err_u):.4e}")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, figsize=(7,7), sharex=True)
        axs[0].plot(x, h_ex, 'k-', label='Ritter (exact, frictionless)')
        axs[0].plot(x, h_num, 'b.', markersize=2, label='numerical')
        axs[0].set_ylabel('h (m)')
        axs[0].legend()
        axs[0].set_title(f'Dam break at t={t:.3f}s')
        axs[1].plot(x, u_ex, 'k-', label='Ritter (exact)')
        axs[1].plot(x, u_num, 'b.', markersize=2, label='numerical')
        axs[1].set_ylabel('u (m/s)')
        axs[1].set_xlabel('x (m)')
        axs[1].legend()
        fig.tight_layout()
        fig.savefig('dam_break_vs_ritter.png', dpi=120)
        print("Saved dam_break_vs_ritter.png")
    except ImportError:
        print("matplotlib not available - skipping plot")


if __name__ == '__main__':
    import sys
    main(*sys.argv[1:])
