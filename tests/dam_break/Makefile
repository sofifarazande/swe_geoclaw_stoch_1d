# Makefile for the dam-break validation test.
# Shared sources are symlinks to the parent directory (see
# tests/lake_at_rest/Makefile).

CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = geoclaw
EXE = xgeo
SETRUN_FILE = setrun.py
OUTDIR = _output
SETPLOT_FILE = setplot.py
PLOTDIR = _plots

OVERWRITE ?= True
RESTART ?= False
FFLAGS ?= -O2

FROM_CLAW := from_claw

COMMON_MODULES = \
  $(FROM_CLAW)/claw_module.f90 \
  $(FROM_CLAW)/utility_module.f90 \
  $(FROM_CLAW)/geoclaw_module.f90 \
  $(FROM_CLAW)/gauges_module.f90 \
  $(FROM_CLAW)/grid_module.f90 \
  $(FROM_CLAW)/topo_module.f90 \

COMMON_SOURCES = \
  $(FROM_CLAW)/claw1.f \
  $(FROM_CLAW)/claw1ez.f \
  $(FROM_CLAW)/setaux.f90 \
  $(FROM_CLAW)/driver.f90 \
  $(FROM_CLAW)/step1.f90 \
  $(FROM_CLAW)/inlinelimiter.f90 \
  $(FROM_CLAW)/copyq1.f \
  $(FROM_CLAW)/opendatafile.f \

MODULES = \
  exner_params_module.f90 \
  sediment_module.f90 \
  sgn_module.f90 \

SOURCES = \
  rp1.f90 \
  src1.f90 \
  bc1.f90 \
  qinit.f90 \
  setprob.f90 \
  out1.f \

include $(CLAWMAKE)

# Module parameters are compiled into every file that uses them, so
# recompile those files whenever a module changes.
sediment_module.o sgn_module.o: exner_params_module.o
rp1.o src1.o bc1.o qinit.o setprob.o: exner_params_module.o sediment_module.o
src1.o setprob.o: sgn_module.o
