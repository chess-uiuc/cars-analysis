# CARS Analysis Tools

[![MATLAB CI](https://github.com/chess-uiuc/cars-analysis/actions/workflows/matlab-ci.yml/badge.svg)](https://github.com/chess-uiuc/cars-analysis/actions/workflows/matlab-ci.yml)
[![CO2 Dual Pump CI](https://github.com/chess-uiuc/cars-analysis/actions/workflows/co2-ci.yml/badge.svg)](https://github.com/chess-uiuc/cars-analysis/actions/workflows/co2-ci.yml)

MATLAB & Fortran analysis tooling for cars data.

## Layout
- `src/matlab/carsft`  – Matlab source for CARSFT functions
- `src/fortran/co2_2pump` - Fortran-based carsft with CO2 and dual pump
- `src/matlab/seed`  – Temporary: Initial files for testing
- `tests/matlab/carsft` – Tests of CARSFT Matlab codes
- `tests/fortran/co2_2pump` - Tests for the CO2/dual pump code
- `tests/matlab/seed` - Temporary: tests of seed project files

## Quick start (MATLAB)
```matlab
% Add all MATLAB source folders to the path:
addpath(genpath('src/matlab'))

% Run the full MATLAB test suite
runtests('tests/matlab')

% Or just test the seed example:
runtests('tests/matlab/seed')

% Run the seed example by hand:
add_one(41)  % -> 42
```

## Quick start (FORTRAN)
### CO2 / Dual Pump
```bash
make -C src/fortran/co2_2pump -j
tests/fortran/co2_2pump/smoke/run_smoke.sh
```
