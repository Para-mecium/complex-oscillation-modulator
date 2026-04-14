function cacheFile = curve_cache_file(cfg, figureId, featureKind, targetValue)
cacheFile = fullfile(cfg.io.cacheDir, figureId, 'curves', ...
    sprintf('%s_%s.mat', circadian.cache_tag(featureKind), circadian.cache_tag(targetValue)));
end
