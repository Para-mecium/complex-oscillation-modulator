function fig = draw_extended_fig2b()
scriptDir = fileparts(mfilename('fullpath'));
cacheFile = fullfile(scriptDir, 'extended_fig2b_cache.mat');

if ~isfile(cacheFile)
    disp('extended_fig2b_cache.mat missing. Generating deterministic cache...');
    generate_extended_fig2b_data();
end

cache = load(cacheFile);
fig = draw_longevity_modulation_panels(cache);
end
