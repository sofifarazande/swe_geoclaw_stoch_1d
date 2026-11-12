"""
setrun.py for the dam-break test: Ritter's dry-bed dam break (h0 = 0.5 m
behind a dam at x0 = 3 m, flat bed), run to t = 0.5 s so that both the
rarefaction tail and the wetting front stay inside the domain. Sediment
and friction are off, for a direct comparison with Ritter's frictionless
solution (check_dam_break.py).
"""

from make_topo import make_topo, L_channel

mx = 300


def setrun(claw_pkg='geoclaw'):

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw', "Expected claw_pkg = 'geoclaw'"

    num_dim = 1
    rundata = data.ClawRunData(claw_pkg, num_dim)

    make_topo()

    probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
    probdata.add_param('sediment_on', 0)
    probdata.add_param('friction_on', 0)
    probdata.add_param('mu_over_sigma', 0.85)   # unused (sediment_on=0)
    probdata.add_param('cd_dep', 0.06)          # unused (sediment_on=0)
    probdata.add_param('bed_alpha', 0.008)      # unused (sediment_on=0)

    clawdata = rundata.clawdata
    clawdata.num_dim = num_dim
    clawdata.lower[0] = 0.0
    clawdata.upper[0] = L_channel
    clawdata.num_cells[0] = mx

    rundata.grid_data.grid_type = 0
    rundata.grid_data.monitor_fgmax = False
    rundata.grid_data.monitor_runup = False
    rundata.grid_data.monitor_total_zeta = False

    clawdata.num_eqn = 4
    clawdata.num_aux = 2
    clawdata.capa_index = 2

    clawdata.t0 = 0.0
    clawdata.restart = False

    clawdata.output_style = 1
    clawdata.num_output_times = 10
    clawdata.tfinal = 0.5
    clawdata.output_t0 = True

    clawdata.output_format = 'ascii'
    clawdata.output_q_components = 'all'
    clawdata.output_aux_components = 'all'
    clawdata.output_aux_onlyonce = True

    clawdata.verbosity = 0

    clawdata.dt_variable = True
    clawdata.dt_initial = 1.e-4
    clawdata.dt_max = 1.0
    clawdata.cfl_desired = 0.45
    clawdata.cfl_max = 0.9
    clawdata.steps_max = 100000

    clawdata.order = 2
    clawdata.num_waves = 2
    clawdata.limiter = [4, 4]
    clawdata.use_fwaves = True
    clawdata.source_split = 'godunov'

    clawdata.num_ghost = 2
    clawdata.bc_lower[0] = 'wall'
    clawdata.bc_upper[0] = 'wall'

    topo_data = rundata.topo_data
    topo_data.topofiles.append([1, 'flat.topo'])

    geo_data = rundata.geo_data
    geo_data.dry_tolerance = 1.e-6
    geo_data.friction_forcing = False
    geo_data.coordinate_system = 1

    rundata.gaugedata.gauges = []

    return rundata


if __name__ == '__main__':
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
