# Coupled SWE / Exner / stochastic-sediment solver on GeoClaw 1D classic.
#
# from_claw/ holds copies of the GeoClaw/Clawpack sources we use, with
# unmodified code (some long comment blocks shortened). Copies, so that
# build products are never written into $CLAW.
# Local files replace GeoClaw's rp1, src1, bc1, b4step1, qinit, setprob, out1.

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

#-------------------------------------------------------------------
include $(CLAWMAKE)

# Module parameters are compiled into every file that uses them, so
# recompile those files whenever a module changes.
sediment_module.o sgn_module.o: exner_params_module.o
rp1.o src1.o bc1.o qinit.o setprob.o: exner_params_module.o sediment_module.o
src1.o setprob.o: sgn_module.o

# Build the rigid-floor topofile z_base.topo
.PHONY: topo
topo:
	$(CLAW_PYTHON) make_topo.py

# Origin of each from_claw/ file under $(CLAW)
GEO_1D  := $(CLAW)/geoclaw/src/1d_classic/shallow
CLASSIC := $(CLAW)/classic/src/1d
FROM_GEO     := claw_module.f90 utility_module.f90 geoclaw_module.f90 gauges_module.f90 \
                grid_module.f90 topo_module.f90 claw1.f claw1ez.f setaux.f90
FROM_CLASSIC := driver.f90 step1.f90 inlinelimiter.f90 copyq1.f opendatafile.f

# Check that the code in from_claw/ matches $(CLAW). Some long comment
# blocks were shortened, so comment and blank lines are ignored
# (lines starting with c/C/* in fixed-form .f files, ! in any file).
STRIP = sed -E -e '/^[[:space:]]*$$/d' -e '/^[[:space:]]*!/d'
.PHONY: verify-from-claw
verify-from-claw:
	@for f in $(FROM_GEO) $(FROM_CLASSIC); do \
	  case " $(FROM_GEO) " in *" $$f "*) src=$(GEO_1D)/$$f ;; *) src=$(CLASSIC)/$$f ;; esac; \
	  case $$f in *.f) extra="-e /^[cC*]/d" ;; *) extra="" ;; esac; \
	  $(STRIP) $$extra $(FROM_CLAW)/$$f > .vfc_a; $(STRIP) $$extra $$src > .vfc_b; \
	  cmp -s .vfc_a .vfc_b || { echo "code differs: $$f"; rm -f .vfc_a .vfc_b; exit 1; }; \
	done; rm -f .vfc_a .vfc_b
	@echo "from_claw/ code matches \$$CLAW (comments ignored)."

# Re-copy from_claw/ after a Clawpack upgrade (restores the full stock comments)
.PHONY: update-from-claw
update-from-claw:
	for f in $(FROM_GEO);     do cp $(GEO_1D)/$$f  $(FROM_CLAW)/; done
	for f in $(FROM_CLASSIC); do cp $(CLASSIC)/$$f $(FROM_CLAW)/; done
