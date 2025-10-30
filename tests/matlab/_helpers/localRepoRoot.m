function root = localRepoRoot()
    here = fileparts(mfilename('fullpath'));
    d = here;
    for i = 1:10
        if isfolder(fullfile(d,'src')) && isfolder(fullfile(d,'tests'))
            root = d; return
        end
        parent = fileparts(d);
        if strcmp(parent, d), break; end
        d = parent;
    end
    ws = getenv('GITHUB_WORKSPACE');
    if ~isempty(ws) && isfolder(ws), root = ws; return; end
    root = pwd;
end
