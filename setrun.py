"""
setrun.py -- main case: a dry, sloped channel with an initial sediment
layer, fed by a constant inflow discharge and sediment flux upstream.

q = (h, hu, z, b)  ->  num_eqn = 4
2-wave HLL f-wave Riemann solver (rp1.f90)  ->  num_waves = 2
aux = (z_base, capacity)  ->  num_aux = 2

L_channel and S0 come from make_topo.py and must match
exner_params_module.f90.
"""

import os
from make_topo import make_topo, L_channel, S0

mx = 600   # number of grid cells


#------------------------------
def setrun(claw_pkg='geoclaw', mu_over_sigma=0.85, cd_dep=0.06, bed_alpha=0.0,
           mann_n_flow=0.025, sgn_on=0, sgn_bed_sigma=0.0):
#------------------------------
    """
    mu_over_sigma, cd_dep, bed_alpha, mann_n_flow, sgn_on, sgn_bed_sigma
    are written to setprob.data and read into the same-named Fortran
    variables; the sweep scripts override them. mann_n_flow is the
    Manning n of the flow friction only (the sediment closure keeps
    mann_n); sgn_on adds the SGN non-hydrostatic correction, whose bed
    is smoothed by a Gaussian of std sgn_bed_sigma [m] (0 = unsmoothed).
    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw', "Expected claw_pkg = 'geoclaw'"

    num_dim = 1
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # Rigid-floor topofile z_base.topo
    make_topo()

    #------------------------------------------------------------------
    # Problem-specific parameters, written to setprob.data and read by
    # setprob.f90 in this order.
    #------------------------------------------------------------------
    probdata = rundata.new_UserData(name='probdata', fname='setprob.data')
    probdata.add_param('sediment_on', 1)
    probdata.add_param('friction_on', 1)
    probdata.add_param('mu_over_sigma', mu_over_sigma)
    probdata.add_param('cd_dep', cd_dep)
    probdata.add_param('bed_alpha', bed_alpha)
    probdata.add_param('mann_n_flow', mann_n_flow)
    probdata.add_param('sgn_on', sgn_on)
    probdata.add_param('sgn_bed_sigma', sgn_bed_sigma)

    #------------------------------------------------------------------
    # Standard Clawpack parameters, written to claw.data:
    #------------------------------------------------------------------
    clawdata = rundata.clawdata

    clawdata.num_dim = num_dim
    clawdata.lower[0] = 0.0
    clawdata.upper[0] = L_channel
    clawdata.num_cells[0] = mx

    # Uniform grid
    rundata.grid_data.grid_type = 0
    rundata.grid_data.monitor_fgmax = False
    rundata.grid_data.monitor_runup = False
    rundata.grid_data.monitor_total_zeta = False

    # ---------------
    # Size of system:
    # ---------------
    clawdata.num_eqn = 4    # h, hu, z, b
    clawdata.num_aux = 2    # z_base (static floor), capacity
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------
    clawdata.t0 = 0.0
    clawdata.restart = False

    # -------------
    # Output times:
    #--------------
    clawdata.output_style = 1
    clawdata.num_output_times = 600
    clawdata.tfinal = 300.0            
    clawdata.output_t0 = True

    clawdata.output_format = 'ascii'
    clawdata.output_q_components = 'all'
    clawdata.output_aux_components = 'all'
    clawdata.output_aux_onlyonce = True

    clawdata.verbosity = 0

    # --------------
    # Time stepping:
    # --------------
    clawdata.dt_variable = True
    clawdata.dt_initial = 1.e-3
    clawdata.dt_max = 1.0

    # Desired / maximum Courant number
    clawdata.cfl_desired = 0.45
    clawdata.cfl_max = 0.9

    clawdata.steps_max = 100000

    # ------------------
    # Method to be used:
    # ------------------
    clawdata.order = 2
    clawdata.num_waves = 2

    # MC limiter on both waves
    clawdata.limiter = [4, 4]

    # rp1 returns f-waves (flux differences), not plain waves
    clawdata.use_fwaves = True

    # Source terms (src1.f90) applied by Godunov splitting
    clawdata.source_split = 'godunov'

    # --------------------
    # Boundary conditions:
    # --------------------
    clawdata.num_ghost = 2

    # 'user' conditions in bc1.f90: inflow on the left, outflow on the right
    clawdata.bc_lower[0] = 'user'
    clawdata.bc_upper[0] = 'user'

    # ---------------
    # Topography:
    # ---------------
    topo_data = rundata.topo_data
    topo_data.topofiles.append([1, 'z_base.topo'])

    # Friction is applied in src1.f90, so GeoClaw's own friction is off
    geo_data = rundata.geo_data
    geo_data.dry_tolerance = 1.e-6
    geo_data.friction_forcing = False
    geo_data.coordinate_system = 1

    # ---------------
    # Gauges: none
    # ---------------
    rundata.gaugedata.gauges = []

    return rundata


if __name__ == '__main__':
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
