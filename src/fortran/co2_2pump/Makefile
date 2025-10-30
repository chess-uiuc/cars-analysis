# === Makefile for co2_2pmp ===
# Dir layout (relative to this Makefile):
#   src/   -> *.for sources
#   data/  -> package-specific inputs (cars.mol, co2.mol, cars.par)
#   build/ -> object files (generated)
#   bin/   -> executables (generated)
MAKEFLAGS += -rR
.SUFFIXES:        # clear implicit rules

# Package name (used to name the binary)
PKG     ?= co2

# Toolchain (override on CLI if needed)
# Compiler: auto-detect, but allow user override (env or CLI)
override FC = gfortran
F77     ?= $(FC)
FFLAGS  = -O2 -ffixed-line-length-none
# Static runtime libs for portability (fine on Linux; macOS ignores static system libs)
# LDFLAGS = -static-libquadmath -static-libgfortran -static-libgcc
# Base static runtime libs
LDFLAGS = -static-libgfortran -static-libgcc

# Probe if the compiler accepts -static-libquadmath; add it if supported.
# (Compiles a tiny program from stdin; no files created.)
QM_FLAG_SUPPORTED := $(shell printf 'program p\nend\n' | \
	$(FC) -x f95 - -o /dev/null -static-libquadmath >/dev/null 2>&1 && echo yes || echo no)

ifeq ($(QM_FLAG_SUPPORTED),yes)
  LDFLAGS += -static-libquadmath
endif
# else
#  # Fallback: still link quadmath if needed (dynamic or static if available)
#  LDFLAGS += -lquadmath
# endif


# If BLAS/LAPACK needed later:
# LDFLAGS += -llapack -lblas

# Directories
SRC_DIR := src
DATA_DIR:= data
BUILD   := build
BIN     := bin
EXE     := $(BIN)/carsfit_$(PKG)

# Sources/objects
SOURCES := $(wildcard $(SRC_DIR)/*.f) \
           $(wildcard $(SRC_DIR)/*.for) \
           $(wildcard $(SRC_DIR)/*.f90)

OBJECTS := $(patsubst $(SRC_DIR)/%,$(BUILD)/%,$(SOURCES))
OBJECTS := $(OBJECTS:.f=.o)
OBJECTS := $(OBJECTS:.for=.o)
OBJECTS := $(OBJECTS:.f90=.o)

.PHONY: all clean distclean print run

all: $(EXE)

$(EXE): $(OBJECTS) | $(BIN)
	$(FC) $(OBJECTS) -o $@ $(LDFLAGS)

# Compile fixed/free form
$(BUILD)/%.o: $(SRC_DIR)/%.for | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC_DIR)/%.f | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/%.o: $(SRC_DIR)/%.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD) $(BIN):
	mkdir -p $@

print:
	@echo "PKG     = $(PKG)"
	@echo "FC      = $(FC)"
	@echo "FFLAGS  = $(FFLAGS)"
	@echo "LDFLAGS = $(LDFLAGS)"
	@echo "SOURCES = $(SOURCES)"
	@echo "EXE     = $(EXE)"

clean:
	rm -rf $(BUILD)

distclean: clean
	rm -rf $(BIN)

# Convenience: run in-place with data files in current dir
run: $(EXE)
	@echo "Running $(EXE) using local data files..."
	cd $(BIN) && cp ../$(DATA_DIR)/* . && ./$(notdir $(EXE))
