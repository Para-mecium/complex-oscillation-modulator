function data = cache_get_or_create(cacheFile, builder)
if isfile(cacheFile)
    loaded = load(cacheFile);
    if isfield(loaded, 'data')
        data = loaded.data;
    else
        data = loaded;
    end
    return
end

cacheDir = fileparts(cacheFile);
if ~isfolder(cacheDir)
    mkdir(cacheDir);
end

data = builder();
save(cacheFile, 'data', '-v7.3');
end
