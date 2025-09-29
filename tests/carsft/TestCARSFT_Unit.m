classdef TestCARSFT_Unit < matlab.unittest.TestCase
    % Focus: small, deterministic checks that don't depend on DB contents

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
        function assignTemps_equilibrium(test)
            [Tr, Tv] = assign_temps(300, 3);
            test.verifyEqual(Tr, 300);
            test.verifyEqual(Tv, [300; 300; 300]);
        end

        function assignTemps_nonEq_defaults_to_Tr(test)
            [Tr, Tv] = assign_temps([400 800], 3);
            test.verifyEqual(Tr, 400);
            test.verifyEqual(Tv, [800; 400; 400]);
        end

        function voigt_matches_lorentz_in_tiny_doppler_limit(test)
            chiamp = 1.0;
            w0     = 2300;
            w      = linspace(2297, 2303, 1201);
            gamma  = 0.2;
            T      = 10;
            cmass  = 1e9;    % absurd mass -> Doppler ~ 0

            chiV = voigt_model(chiamp, w0, w, gamma, T, cmass);
            chiL = chiamp ./ ((w0 - w) - 1i*gamma/2);

            chiV = chiV ./ max(abs(chiV));
            chiL = chiL ./ max(abs(chiL));

            rel = norm(chiV - chiL) / norm(chiL);
            test.verifyLessThan(rel, 0.05);
        end
    end
end
