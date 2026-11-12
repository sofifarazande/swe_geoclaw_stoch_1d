"""
model_params.py - model parameters read directly from the solver, so the
analysis scripts never carry their own copies:

  - Fortran `parameter` constants from exner_params_module.f90 and
    sediment_module.f90 (expressions evaluated in declaration order),
  - run-time values from setprob.data (mu_over_sigma, cd_dep, bed_alpha, ...),
  - grid spacing dx from claw.data.

Usage:
    from model_params import load
    P = load()            # P['Q_in'], P['d50'], P['mu_over_sigma'], P['dx'], ...
    P = load(outdir)      # .data values of that run, if saved in outdir
    python scripts/model_params.py [OUTDIR]     # print all values
"""
import math
import os
import re

CASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FORTRAN_MODULES = ('exner_params_module.f90', 'sediment_module.f90')

_DECL = re.compile(r'^\s*(real|integer|logical)\b[^:]*\bparameter\b[^:]*::\s*(\w+)\s*=\s*(.+)$',
                   re.IGNORECASE)


def _to_python(expr):
    expr = expr.split('!')[0].strip()
    expr = re.sub(r'(\d)[dD]([+-]?\d)', r'\1e\2', expr)        # 1.0d-30 -> 1.0e-30
    expr = re.sub(r'\.true\.', 'True', expr, flags=re.IGNORECASE)
    expr = re.sub(r'\.false\.', 'False', expr, flags=re.IGNORECASE)
    return expr


def read_fortran_parameters(paths):
    """Evaluate every `parameter` declaration in the given Fortran files."""
    ns = {'sqrt': math.sqrt, 'log': math.log, 'log10': math.log10, 'exp': math.exp,
          'sin': math.sin, 'cos': math.cos, 'atan': math.atan, 'abs': abs,
          'max': max, 'min': min, 'nint': round}
    params = {}
    for path in paths:
        with open(path) as f:
            for line in f:
                m = _DECL.match(line)
                if not m:
                    continue
                name, expr = m.group(2), _to_python(m.group(3))
                value = eval(expr, ns)
                if m.group(1).lower() == 'integer':
                    value = int(value)
                ns[name] = value
                params[name] = value
    return params


def read_clawpack_data(path):
    """Values of a Clawpack .data file ('value  =: name' lines)."""
    values = {}
    with open(path) as f:
        for line in f:
            if '=:' not in line:
                continue
            value, name = (s.strip() for s in line.split('=:', 1))
            try:
                values[name] = int(value)
            except ValueError:
                try:
                    values[name] = float(value)
                except ValueError:
                    values[name] = value
    return values


def load(outdir=None, case_dir=CASE_DIR):
    """Fortran parameters + setprob.data + dx (from claw.data). The .data
    files are taken from outdir when it has them (the values that run
    used), else from case_dir (the values the next run will use)."""
    P = read_fortran_parameters([os.path.join(case_dir, f) for f in FORTRAN_MODULES])
    data_dir = case_dir
    if outdir and os.path.exists(os.path.join(outdir, 'setprob.data')):
        data_dir = outdir
    P.update(read_clawpack_data(os.path.join(data_dir, 'setprob.data')))
    claw = read_clawpack_data(os.path.join(data_dir, 'claw.data'))
    P['dx'] = (claw['upper'] - claw['lower']) / claw['num_cells']
    return P


def settling_velocity(P):
    """Dietrich (1982) settling velocity, as sediment_module.f90's ws_settle."""
    nu_k = 1.0e-6
    d_star = P['s_rel'] * P['g'] * P['d50']**3 / nu_k**2
    L = math.log10(d_star)
    logF = -3.76715 + 1.92944*L - 0.09815*L**2 - 0.00575*L**3 + 0.00056*L**4
    return (P['s_rel'] * P['g'] * nu_k)**(1.0/3.0) * (10.0**logF)**(1.0/3.0)


def normal_flow(P):
    """Uniform (normal) flow for the inlet discharge: Manning friction slope
    n^2 u^2 / h^(4/3) = S0 with u*h = Q_in/B_channel. Returns (h, u)."""
    q = P['Q_in'] / P['B_channel']
    h = (P['mann_n'] * q / math.sqrt(P['S0']))**0.6
    return h, q / h


def flow_from_output(outdir, dry_tol=1.e-6):
    """Mean wet depth and velocity over the second half of a run. Returns (h, u)."""
    import numpy as np
    from plot_results import available_frames, read_frame
    frames = available_frames(outdir)
    hs, us = [], []
    for n in frames[len(frames)//2:]:
        x, h, hu, z, b, eta, t = read_frame(outdir, n)
        wet = h > dry_tol
        if np.any(wet):
            hs.append(h[wet].mean())
            us.append((hu[wet] / h[wet]).mean())
    return float(np.mean(hs)), float(np.mean(us))


def base_flow(P, h0=None, u0=None, from_output=None):
    """(h0, u0) for the linear analyses: explicit values if given, else the
    mean late-time flow of a run (from_output), else the normal flow."""
    if h0 is not None and u0 is not None:
        return h0, u0, 'given'
    if from_output:
        return (*flow_from_output(from_output), f'mean of 2nd half of {from_output}')
    return (*normal_flow(P), 'normal flow for Q_in')


if __name__ == '__main__':
    import sys
    P = load(sys.argv[1] if len(sys.argv) > 1 else None)
    for k, v in P.items():
        print(f'{k:>18} = {v}')
    h, u = normal_flow(P)
    print(f'\nnormal flow: h = {h:.5f} m, u = {u:.4f} m/s, '
          f'Fr = {u/math.sqrt(P["g"]*h):.3f};  ws = {settling_velocity(P)*1000:.2f} mm/s')
