function ensure_paths()
persistent isReady
if ~isempty(isReady) && isReady
    return
end

packageDir = fileparts(mfilename('fullpath'));
circadianDir = fileparts(packageDir);
fmamDir = fileparts(circadianDir);
matcontDir = fullfile(fmamDir, 'MatCont7p6');

addpath(fmamDir);
addpath(circadianDir);
addpath(fullfile(fmamDir, 'PO_extract'));
addpath(genpath(matcontDir));

isReady = true;
end
