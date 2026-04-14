function cacheFile = curve_cache_file(cfg, figureId, featureKind, targetValue)
cacheFile = fullfile(cfg.io.cacheDir, figureId, 'curves', ...
    sprintf('%s_%s.mat', flexmod.cache_tag(featureKind), flexmod.cache_tag(targetValue)));
end
