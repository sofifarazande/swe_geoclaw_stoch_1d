"""
Quantitative well-balancing check: scans every fort.qNNNN output frame
and reports max|u| and max|eta-eta0| over the whole domain, as a
function of time. Run after `make .output` (or `make .plots`).
"""
import glob
import re
import numpy as np

def read_frame(fname):
    data = np.loadtxt(fname, skiprows=5)
    h, hu, z = data[:,0], data[:,1], data[:,2]
    return h, hu, z

def read_time(fname):
    tname = fname.replace('fort.q', 'fort.t')
    with open(tname) as f:
        return float(f.readline().split()[0])

def main(outdir='_output'):
    qfiles = sorted(glob.glob(f'{outdir}/fort.q[0-9][0-9][0-9][0-9]'))
    if not qfiles:
        print(f'No fort.q files found in {outdir} - run the case first (make .output).')
        return

    h0, hu0, z0 = read_frame(qfiles[0])
    eta0 = h0 + z0

    print(f"{'t':>10} {'max|u|':>14} {'max|eta-eta0|':>16}")
    worst_u = 0.0
    worst_eta = 0.0
    for fname in qfiles:
        h, hu, z = read_frame(fname)
        t = read_time(fname)
        eta = h + z
        wet = h > 1.e-6
        u = np.zeros_like(h)
        u[wet] = hu[wet] / h[wet]
        max_u = np.max(np.abs(u))
        max_deta = np.max(np.abs(eta - eta0))
        worst_u = max(worst_u, max_u)
        worst_eta = max(worst_eta, max_deta)
        print(f"{t:10.3f} {max_u:14.6e} {max_deta:16.6e}")

    print()
    print(f"Worst max|u|  over the whole run: {worst_u:.6e}")
    print(f"Worst max|eta-eta0| over the whole run: {worst_eta:.6e}")
    print()
    if worst_u < 1.e-10:
        print("=> Essentially machine-precision: the scheme is well-balanced.")
    elif worst_u < 1.e-4:
        print("=> Small but not machine-precision: mild imbalance from the")
        print("   point-wise bed-slope source; likely fine in practice but a")
        print("   hydrostatic-reconstruction upgrade would remove it.")
    else:
        print("=> Not well-balanced: velocities are growing to a physically")
        print("   significant size at rest. The point-wise bed-slope source")
        print("   needs the hydrostatic-reconstruction upgrade.")

if __name__ == '__main__':
    import sys
    main(*sys.argv[1:])
