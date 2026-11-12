"""Bump-bed topofile for the lake-at-rest test."""

L_channel = 6.0
eta_const = 0.35   # flat free-surface level imposed by qinit.f90 (must stay > bump peak)

def bump(x):
    return 0.1 + 0.2 * 2.718281828459045 ** (-((x - 3.0) / 0.6) ** 2)

def make_topo(fname='bump.topo', npts=201):
    with open(fname, 'w') as f:
        f.write('%d\n' % npts)
        for i in range(npts):
            x = L_channel * i / (npts - 1)
            f.write('%.10f %.10f\n' % (x, bump(x)))

if __name__ == '__main__':
    make_topo()
