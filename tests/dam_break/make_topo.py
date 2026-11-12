"""Flat bed (z_base=0 everywhere) for the dam-break test."""

L_channel = 6.0

def make_topo(fname='flat.topo'):
    with open(fname, 'w') as f:
        f.write('2\n')
        f.write('%.10f %.10f\n' % (0.0, 0.0))
        f.write('%.10f %.10f\n' % (L_channel, 0.0))

if __name__ == '__main__':
    make_topo()
