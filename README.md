# CARS Analysis Tools

[![MATLAB CI](https://github.com/chess-uiuc/cars-analysis/actions/workflows/matlab-ci.yml/badge.svg)](https://github.com/chess-uiuc/cars-analysis/actions/workflows/matlab-ci.yml)
[![CO2 Dual Pump CI](https://github.com/chess-uiuc/cars-analysis/actions/workflows/co2-ci.yml/badge.svg)](https://github.com/chess-uiuc/cars-analysis/actions/workflows/co2-ci.yml)

MATLAB & Fortran analysis tooling for cars data.

## Layout

### Main tools
- `src/matlab/carsft`  – Matlab source for CARSFT functions
- `src/matlab/+cfx` - Matlab interface for Fortran carsfit (CO2/2pump)
- `src/fortran/co2_2pump` - Fortran-based carsfit (CO2/2pump)
- `src/fortran/co2_2pump/data` - Data files for the fortran-based exe
- `src/python/carsfit-tools` - Python utilities for carsfit (e.g. plotting)

### Tests
- `tests/matlab/carsft` – Tests of CARSFT Matlab codes
- `tests/matlab/cfx` - Tests of the interface to Fortran carsfit (CO2/2pump)
- `tests/fortran/co2_2pump` - Tests for the CO2/dual pump code

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

## Quick start (MATLAB/FORTRAN carsfit Interface)

1) Build or obtain the binary (See FORTRAN below):
   `src/fortran/co2_2pump/bin/carsfit_co2`

2) In MATLAB:

```matlab
addpath(genpath('src/matlab'));

% Run carsfit_co2 from MATLAB (8 Enters, then N, N)
seq = [repmat("",8,1); "N"; "N"];
out = cfx.run_carsfit_script(seq, struct('workdir',"runs/demo1"));

% Load and plot the primary CSV
T = out.tables(out.primary_csv);
figure; plot(T{:,1}, T{:,2}, 'LineWidth', 1.2); grid on; xlabel('X'); ylabel('Y');

% Optional: generate a PNG plot via the Python tool
R = cfx.plot_csv_with_python(fullfile(out.workdir, out.primary_csv));
disp("PNG saved as " + R.png_out);
```

## Quick start (FORTRAN)
### CO2 / Dual Pump
```bash
make -C src/fortran/co2_2pump -j
tests/fortran/co2_2pump/smoke/run_smoke.sh
```
