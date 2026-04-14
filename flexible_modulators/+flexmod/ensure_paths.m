function ensure_paths()
persistent isReady
if ~isempty(isReady) && isReady
    return
end

packageDir = fileparts(mfilename('fullpath'));
flexDir = fileparts(packageDir);
fmamDir = fileparts(flexDir);
matcontDir = fullfile(fmamDir, 'MatCont7p6');

addpath(fmamDir);
addpath(flexDir);
addpath(fullfile(fmamDir, 'PO_extract'));
addpath(genpath(matcontDir));

isReady = true;
end
