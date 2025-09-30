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
        function scalar_T_equals_vector_equilibrium_T(test)
            % Passing T as a scalar vs a length-1 vector should be equivalent.
            wexp  = [2200:0.1:2400];
            X     = [1 0 0 0];
            P     = 1.0; dtp = 5.0; dtau3 = 0.0; alpha = 0.0; dwp = 1.0;

            [S1, chi1, w1] = CARSFT_dev(wexp, 300,  P, X, dtp, dtau3, alpha, dwp);
            [S2, chi2, w2] = CARSFT_dev(wexp, [300], P, X, dtp, dtau3, alpha, dwp);

            test.verifyEqual(w1, w2, "AbsTol", 0);
            test.verifyLessThan(norm(S1 - S2) / max(1,norm(S1)), 1e-12);
            test.verifyLessThan(norm(chi1 - chi2) / max(1,norm(chi1)), 1e-12);
        end
 
        function interpolation_resolution_is_neutral_after_norm(test)
            % Different wexp resolutions should yield the same shape after normalization.
            X=[1 0 0 0]; T=300; P=1.0; dtp=100.0; dtau3=0.0; alpha=0.0; dwp=1.0;

            w_fine = [2000:0.05:2600];
            w_coarse = [2000.0:0.1:2600];

            [Sfine, ~, ~]   = CARSFT_dev(w_fine,   T, P, X, dtp, dtau3, alpha, dwp);
            [Scoarse, ~, ~] = CARSFT_dev(w_coarse, T, P, X, dtp, dtau3, alpha, dwp);

            % Compare on the fine grid by interpolating the coarse result
            Scoarse_up = interp1(w_coarse, Scoarse, w_fine, 'linear', 'extrap');

            fprintf('Sfine:   min=%g max=%g mean=%g NaN? %d Inf? %d\n', ...
                min(Sfine), max(Sfine), mean(Sfine), any(isnan(Sfine)), any(isinf(Sfine)));

            fprintf('Scoarse: min=%g max=%g mean=%g NaN? %d Inf? %d\n', ...
                min(Scoarse), max(Scoarse), mean(Scoarse), any(isnan(Scoarse)), any(isinf(Scoarse)));

            fprintf('Scoarse_up: min=%g max=%g mean=%g NaN? %d Inf? %d\n', ...
                min(Scoarse_up), max(Scoarse_up), mean(Scoarse_up), ...
                any(isnan(Scoarse_up)), any(isinf(Scoarse_up)));

            % Both are normalized; allow small interp error
            rel = norm(Sfine - Scoarse_up) / max(1, norm(Sfine));
            test.verifyLessThan(rel, 0.02);  % 2% tolerance
        end
    end
end
