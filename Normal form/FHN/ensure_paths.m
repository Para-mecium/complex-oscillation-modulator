function ensure_paths()
persistent isReady
if ~isempty(isReady) && isReady
    return
end

fhnDir = fileparts(mfilename('fullpath'));
normalFormDir = fileparts(fhnDir);
fmamDir = fileparts(normalFormDir);

addpath(fmamDir);
addpath(normalFormDir);
addpath(fhnDir);

isReady = true;
end
