classdef TestCarsftRepoSmoke < matlab.unittest.TestCase
    % Minimal, structure-safe smoke tests for CARSFT.
    % These don't assume any specific function names.

    properties (Access = private)
      CarsftDir % instead now set below 
    end
    methods (TestClassSetup)
        function addRepoPaths(tc)
            % 1) Make the shared helper visible
            here     = fileparts(mfilename('fullpath'));   % .../tests/carsft
            helpers  = fullfile(here, '..', '_helpers');   % .../tests/_helpers
            addpath(helpers, '-begin');

            % 2) Resolve repo root using the centralized helper
            repoRoot     = localRepoRoot();
            tc.CarsftDir = fullfile(repoRoot, 'src', 'matlab', 'carsft');

            % 3) Put src/carsft on the MATLAB path for the tests
            addpath(genpath(tc.CarsftDir));
        end
    end

    methods (Test)
        function testCarsftFolderExistsAndHasCode(tc)
            tc.assertTrue(isfolder(tc.CarsftDir), ...
                'Missing src/carsft directory.');
            mfiles = dir(fullfile(tc.CarsftDir, '**', '*.m'));
            tc.assertGreaterThanOrEqual(numel(mfiles), 1, ...
                'No .m files found in src/carsft.');
        end

        function testCarsftFilesAreOnPath(tc)
            % Pick a few top-level function files and make sure MATLAB can resolve them.
            mfiles = dir(fullfile(tc.CarsftDir, '*.m'));
            if isempty(mfiles)
                tc.assumeFail('No top-level .m files in src/carsft to resolve.');
            end
            % Check up to 5 files to keep the test quick.
            for k = 1:min(5, numel(mfiles))
                [~, fname] = fileparts(mfiles(k).name);
                % Skip private/package/class folders by construction (top-level only).
                resolved = which(fname);
                tc.assertNotEmpty(resolved, sprintf('Function "%s" is not on path.', fname));
            end
        end

        function testZeroArgEntrypointsRun(tc)
            % Call any zero-argument top-level functions (if any) and
            % assert they don't error. This is a pure smoke test.
            mfiles = dir(fullfile(tc.CarsftDir, '*.m'));
            ranAny = false;

            for k = 1:numel(mfiles)
                [~, fname] = fileparts(mfiles(k).name);

                % Skip scripts: we only want functions (best-effort check)
                % Heuristic: read first non-empty, non-comment line and see if it starts with "function"
                isFunction = localIsFunctionFile(fullfile(mfiles(k).folder, mfiles(k).name));
                if ~isFunction
                    continue
                end

                fh = str2func(fname);
                % nargin for function handle works for most simple functions
                try
                    nin = nargin(fh);
                catch
                    % If MATLAB can't determine nargin, skip conservatively.
                    continue
                end
                if nin == 0
                    ranAny = true;
                    tc.verifyWarningFree(@() fh(), ...
                        sprintf('Zero-arg function "%s" threw or warned.', fname));
                end
            end

            % Not all libraries have zero-arg entrypoints; don't fail the suite for that,
            % just mark the assumption so CI shows it as "skipped".
            tc.assumeTrue(ranAny, ...
                'No zero-argument top-level functions found to smoke-test.');
        end
    end
end

function tf = localIsFunctionFile(fpath)
% Return true if file appears to be a function (not a script).
    tf = false;
    try
        txt = fileread(fpath);
    catch
        return
    end
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    for i = 1:numel(lines)
        L = strtrim(lines{i});
        if isempty(L) || startsWith(L,'%')
            continue
        end
        tf = startsWith(L, 'function');
        return
    end
end
