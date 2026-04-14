function fig = draw_figS1a()
scriptDir = fileparts(mfilename('fullpath'));
cacheFile = fullfile(scriptDir, 'figS1a_cache.mat');

if ~isfile(cacheFile)
    disp('figS1a_cache.mat missing. Generating deterministic cache...');
    generate_figS1a_data();
end

cache = load(cacheFile);
fig = draw_longevity_modulation_panels(cache);
end
