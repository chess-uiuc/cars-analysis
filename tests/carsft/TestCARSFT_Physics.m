classdef TestCARSFT_Physics < matlab.unittest.TestCase
    % Focus: real pipeline with real CARS.mol; physics-based assertions

    methods (TestClassSetup)
          function addRepoPaths(tc)
             % 1. Add the _helpers folder so localRepoRoot is visible
             here = fileparts(mfilename('fullpath'));             % e.g. tests/carsft
             helpers = fullfile(here, '..', '_helpers');
             addpath(helpers, '-begin');

             % 2. Add the carsft source directory from the repo root
             repoRoot = localRepoRoot();
             addpath(genpath(fullfile(repoRoot, 'src', 'carsft')));
           end
    end

    methods (Test)
        function e2e_runs_and_has_structure(test)
            wexp  = [2050:0.1:2350];
            T     = 3000;  P = 0.2;  X = [1 0.5 0 0];
            dtp   = 8000;  dtau3 = 0.0;  alpha = 0.0;  dwp = 1.0;
            [S, chi, w] = CARSFT_dev(wexp, T, P, X, dtp, dtau3, alpha, dwp);

            test.verifySize(S, size(wexp));
            test.verifyFalse(any(~isfinite(S)));
            test.verifyGreaterThanOrEqual(min(S), 0);
            test.verifyEqual(max(S), 1, "AbsTol", 1e-12);
            test.verifySize(chi, [1 numel(w)]);
            test.verifyGreaterThan(std(S), 1e-6); % not flat
        end

        function chi_scales_with_pressure_over_temperature(test)
            wexp  = [2050:0.1:2350];
            X     = [1 0 0 0];
            dtp   = 1000; dtau3 = 0.0; alpha = 0.0; dwp = 1.0;

            T1=300; P1=1.0; [~, chi1, ~] = CARSFT_dev(wexp, T1, P1, X, dtp, dtau3, alpha, dwp);
            T2=300; P2=2.0; [~, chi2, ~] = CARSFT_dev(wexp, T2, P2, X, dtp, dtau3, alpha, dwp);
            r = max(abs(chi2)) / max(abs(chi1));
            test.verifyEqual(r, 2.0, "RelTol", 0.10);

            T3=600; P3=1.0; [~, chi3, ~] = CARSFT_dev(wexp, T3, P3, X, dtp, dtau3, alpha, dwp);
            rT = max(abs(chi3)) / max(abs(chi1));
            test.verifyEqual(rT, 0.5, "RelTol", 0.10);
        end

        function probe_and_instrument_broaden_lines(test)
            wexp  = [2200:0.1:2400];
            T=300; P=1.0; X=[1 0 0 0]; dtau3=0.0; alpha=0.0;

            dtp=1000.0; dwp1=.5; dwp2=1.0;
            [S1, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp, dtau3, alpha, dwp1);
            [S2, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp, dtau3, alpha, dwp2);

            f1 = fwhm(wexp, S1);
            f2 = fwhm(wexp, S2);
            test.verifyGreaterThan(f2, f1);

            dwp=0.5; dtp_long=1000.0; dtp_short=500.0;
            [S3, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp_long, dtau3, alpha, dwp);
            [S4, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp_short, dtau3, alpha, dwp);

            f3 = fwhm(wexp, S3);
            f4 = fwhm(wexp, S4);
            test.verifyGreaterThan(f4, f3);
        end
    end
end

% --- helper: crude FWHM on normalized S ---
function width = fwhm(x, y)
    y = y(:); x = x(:);
    [ymax, imax] = max(y);
    if ~isfinite(ymax) || ymax <= 0
        width = NaN; return;
    end
    half = ymax/2;
    iL = find(y(1:imax) <= half, 1, 'last');
    if isempty(iL), xL = x(1); else
        xL = interp1(y(iL:iL+1), x(iL:iL+1), half, 'linear', 'extrap');
    end
    iR = imax - 1 + find(y(imax:end) <= half, 1, 'first');
    if isempty(iR), xR = x(end); else
        xR = interp1(y(iR-1:iR), x(iR-1:iR), half, 'linear', 'extrap');
    end
    width = xR - xL;
end
