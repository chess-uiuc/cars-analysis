classdef TestCarsftSmoke < matlab.unittest.TestCase
    % Minimal, structure-safe smoke tests for CARSFT.
    % These don't assume any specific function names.

    properties (Constant)
        CarsftDir = fullfile(localRepoRoot(), 'src', 'carsft');
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

function root = localRepoRoot()
% Find the repository root robustly in CI and locally.
    % 1) Try walking up from this test file to find a folder that has src/ and tests/
    here = fileparts(mfilename('fullpath'));
    d = here;
    for i = 1:10   % don't walk forever
        if isfolder(fullfile(d,'src')) && isfolder(fullfile(d,'tests'))
            root = d; return
        end
        parent = fileparts(d);
        if strcmp(parent, d), break; end
        d = parent;
    end
    % 2) GitHub Actions workspace, if present
    ws = getenv('GITHUB_WORKSPACE');
    if ~isempty(ws) && isfolder(ws)
        root = ws; return
    end
    % 3) Fallback: current directory
    root = pwd;
end
