function cacheFile = mark_cache_file(cfg, figureId, markId)
cacheFile = fullfile(cfg.io.cacheDir, figureId, 'marks', ...
    sprintf('%s.mat', flexmod.cache_tag(markId)));
end
