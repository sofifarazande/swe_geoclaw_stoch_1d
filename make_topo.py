"""
Write z_base.topo, the rigid channel floor z_base(x) = S0*(L_channel - x),
which GeoClaw's topo_module reads into aux(1,:). Two points suffice since
z_base is linear.

L_channel and S0 must match exner_params_module.f90.
"""

L_channel = 6.0
S0 = 0.044

def make_topo(fname='z_base.topo'):
    with open(fname, 'w') as f:
        f.write('2\n')
        f.write('%.10f %.10f\n' % (0.0, S0 * L_channel))
        f.write('%.10f %.10f\n' % (L_channel, 0.0))

if __name__ == '__main__':
    make_topo()
