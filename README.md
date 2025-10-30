# chess-cars-analysis

MATLAB-based analysis tooling for CHESS cars data.

![MATLAB CI](https://github.com/MTCam/chess-cars-analysis/actions/workflows/matlab-ci.yml/badge.svg)

## Layout
- `src/matlab/carsft`  – Sean Kearney's original Matlab source
- `src/matlab/seed`  – Temporary: Initial files for testing
- `tests/matlab/carsft` – Tests of CARSFT Matlab codes
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