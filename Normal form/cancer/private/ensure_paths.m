function ensure_paths()
persistent isReady
if ~isempty(isReady) && isReady
    return
end

privateDir = fileparts(mfilename('fullpath'));
cancerDir = fileparts(privateDir);
normalFormDir = fileparts(cancerDir);
fmamDir = fileparts(normalFormDir);

addpath(fmamDir);
addpath(normalFormDir);
addpath(cancerDir);

isReady = true;
end
