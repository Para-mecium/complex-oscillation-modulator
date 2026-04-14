function result = reproduce_all_fig3(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

result = struct();
result.fig3b = reproduce_fig3b(cfg);
result.fig3c = reproduce_fig3c(cfg);
result.fig3d = reproduce_fig3d(cfg);
result.fig3f = reproduce_fig3f(cfg);
end
