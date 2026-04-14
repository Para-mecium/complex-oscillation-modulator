function result = reproduce_all_fig5(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

result = struct();
result.fig5b = reproduce_fig5b(cfg);
result.fig5c = reproduce_fig5c(cfg);
result.fig5d = reproduce_fig5d(cfg);
result.figS2a = reproduce_figS2a(cfg);
end
