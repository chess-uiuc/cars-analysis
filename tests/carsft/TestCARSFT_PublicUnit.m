classdef TestCARSFT_PublicUnit < matlab.unittest.TestCase
    % Public-API tests that exercise invariants via CARSFT_dev only.
    % No direct calls to nested/private helpers required.

    properties (Access = private)
        CarsftDir  % set in TestClassSetup after helpers are on path
    end

    methods (TestClassSetup)
        function addRepoPaths(tc)
            % 1) Make the shared helper visible
            here    = fileparts(mfilename('fullpath'));      % .../tests/carsft
            helpers = fullfile(here, '..', '_helpers');      % .../tests/_helpers
            addpath(helpers, '-begin');

            % 2) Resolve repo root and add src/carsft
            repoRoot     = localRepoRoot();
            tc.CarsftDir = fullfile(repoRoot, 'src', 'carsft');
            addpath(genpath(tc.CarsftDir));
        end
    end

    methods (Test)
        function chi_scales_with_pressure_over_temperature(test)
            % Verify the intended P/T scaling on chi (pre-normalization).
            wexp  = linspace(200, 4000, 4001);
            X     = [1 0 0 0];
            dtp   = 0.2; dtau3 = 0.0; alpha = 0.0; dwp = 0.0;

            T1=300; P1=1.0; [~, chi1, ~] = CARSFT_dev(wexp, T1, P1, X, dtp, dtau3, alpha, dwp);
            T2=300; P2=2.0; [~, chi2, ~] = CARSFT_dev(wexp, T2, P2, X, dtp, dtau3, alpha, dwp);
            rP = max(abs(chi2)) / max(abs(chi1));
            test.verifyEqual(rP, 2.0, "RelTol", 0.10);

            T3=600; P3=1.0; [~, chi3, ~] = CARSFT_dev(wexp, T3, P3, X, dtp, dtau3, alpha, dwp);
            rT = max(abs(chi3)) / max(abs(chi1));
            test.verifyEqual(rT, 0.5, "RelTol", 0.10);
        end

        function scalar_T_equals_vector_equilibrium_T(test)
            % Passing T as a scalar vs a length-1 vector should be equivalent.
            wexp  = linspace(2200, 2400, 2001);
            X     = [1 0 0 0];
            P     = 1.0; dtp = 0.2; dtau3 = 0.0; alpha = 0.0; dwp = 0.0;

            [S1, chi1, w1] = CARSFT_dev(wexp, 300,  P, X, dtp, dtau3, alpha, dwp);
            [S2, chi2, w2] = CARSFT_dev(wexp, [300], P, X, dtp, dtau3, alpha, dwp);

            test.verifyEqual(w1, w2, "AbsTol", 0);
            test.verifyLessThan(norm(S1 - S2) / max(1,norm(S1)), 1e-12);
            test.verifyLessThan(norm(chi1 - chi2) / max(1,norm(chi1)), 1e-12);
        end

        function probe_and_instrument_broaden_lines(test)
            % Shorter probe (smaller dtp) and larger instrument width (dwp)
            % should both broaden the normalized spectrum S.
            wexp  = linspace(2200, 2400, 4001);
            T=300; P=1.0; X=[1 0 0 0]; dtau3=0.0; alpha=0.0;

            % Instrument broadening
            dtp = 2.0; dwp1 = 0.0; dwp2 = 0.8;
            [S1, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp, dtau3, alpha, dwp1);
            [S2, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp, dtau3, alpha, dwp2);
            f1 = fwhm(wexp, S1); f2 = fwhm(wexp, S2);
            test.verifyGreaterThan(f2, f1);

            % Probe broadening
            dwp = 0.0; dtp_long = 2.0; dtp_short = 0.2;
            [S3, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp_long,  dtau3, alpha, dwp);
            [S4, ~, ~] = CARSFT_dev(wexp, T, P, X, dtp_short, dtau3, alpha, dwp);
            f3 = fwhm(wexp, S3); f4 = fwhm(wexp, S4);
            test.verifyGreaterThan(f4, f3);
        end

        function interpolation_resolution_is_neutral_after_norm(test)
            % Different wexp resolutions should yield the same shape after normalization.
            X=[1 0 0 0]; T=300; P=1.0; dtp=0.2; dtau3=0.0; alpha=0.0; dwp=0.2;

            w_fine = linspace(2000, 2600, 8001);
            w_coarse = linspace(2000, 2600, 1601);

            [Sfine, ~, ~]   = CARSFT_dev(w_fine,   T, P, X, dtp, dtau3, alpha, dwp);
            [Scoarse, ~, ~] = CARSFT_dev(w_coarse, T, P, X, dtp, dtau3, alpha, dwp);

            % Compare on the fine grid by interpolating the coarse result
            Scoarse_up = interp1(w_coarse, Scoarse, w_fine, 'linear', 'extrap');

            % Both are normalized; allow small interp error
            rel = norm(Sfine - Scoarse_up) / max(1, norm(Sfine));
            test.verifyLessThan(rel, 0.02);  % 2% tolerance
        end
    end
end

% --- local helper: crude FWHM on normalized S (max(S) ~ 1) ---
function width = fwhm(x, y)
    y = y(:); x = x(:);
    [ymax, imax] = max(y);
    if ~isfinite(ymax) || ymax <= 0
        width = NaN; return;
    end
    half = ymax/2;

    % Left crossing
    iL = find(y(1:imax) <= half, 1, 'last');
    if isempty(iL), xL = x(1); else
        xL = interp1(y(iL:iL+1), x(iL:iL+1), half, 'linear', 'extrap');
    end

    % Right crossing
    iR = imax - 1 + find(y(imax:end) <= half, 1, 'first');
    if isempty(iR), xR = x(end); else
        xR = interp1(y(iR-1:iR), x(iR-1:iR), half, 'linear', 'extrap');
    end

    width = xR - xL;
end
