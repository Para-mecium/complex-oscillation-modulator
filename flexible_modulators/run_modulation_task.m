function result = run_modulation_task(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);
result = flexmod.run_fmam_task(cfg);
end
